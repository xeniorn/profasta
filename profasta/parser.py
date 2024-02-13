"""This module manages parsers for the headers of FASTA records.

This module provides classes for parsing the headers of FASTA records into a structured
format and writing the structured format back to a header string. The default FASTA
header parsers are registered in a global registry, which can be accessed via the
`get_parser` function and the name of the parser. New parsers must be registered via
the `register_parser` function before they become available in the other modules.

Classes:
    AbstractParsedHeader (Protocol): Interface for representing a parsed FASTA header.
    ParsedHeader: Representation of a parsed FASTA header.
    HeaderParser (Protocol): Interface for a FASTA header parser.
    DefaultParser: Default FASTA header parser.
    UniprotParser: Parser for Uniprot FASTA headers.

Functions:
    register_parser: Register a custom FASTA header parser by name.
    get_parser: Get a registered FASTA header parser by name.

Constants:
    PARSER_REGISTRY: Dictionary mapping parser names to header parsers. The built-in
        parsers are registered as "default" and "uniprot" and can be retrieved via the
        `get_parser` function.
"""

from dataclasses import dataclass, field
from typing import Protocol
import re


class AbstractParsedHeader(Protocol):
    """Abstract parsed FASTA header.

    Attributes:
        identifier: The unique identifier of the FASTA entry.
        header: The FASTA header, not containing the starting ">" character.
        header_fields: The parsed header fields as a dictionary.
    """

    identifier: str
    header: str
    header_fields: dict[str, str]


@dataclass
class ParsedHeader:
    """Parsed FASTA header.

    Attributes:
        identifier: The unique identifier of the FASTA entry.
        header: The FASTA header, not containing the starting ">" character.
        header_fields: The parsed header fields as a dictionary.
    """

    identifier: str
    header: str
    header_fields: dict[str, str] = field(default_factory=dict)


class HeaderParser(Protocol):
    """Abstract header parser."""

    @classmethod
    def parse(self, header: str) -> AbstractParsedHeader:
        """Parse a FASTA header string into a ParsedHeader object."""
        ...

    @classmethod
    def write(self, parsed_header: AbstractParsedHeader) -> str:
        """Write a FASTA header string from a ParsedHeader object."""
        ...


class DefaultParser:
    """Default FASTA header parser.

    The `parse` method returns a ParsedHeader object with the identifier being the
    first whitespace-separated word of the header. The rest of the header is stored
    in the "description" field of the `header_fields` dictionary, which might be an
    empty string.

    The `write` method returns the original `header` string from the parsed_header.
    """

    @classmethod
    def parse(cls, header: str) -> ParsedHeader:
        """Parse a FASTA header string into a ParsedHeader object."""
        split_header = header.split(maxsplit=1)
        _id = split_header[0]
        fields = {"description": split_header[1]} if len(split_header) > 1 else {}
        return ParsedHeader(_id, header, fields)

    @classmethod
    def write(cls, parsed_header: AbstractParsedHeader) -> str:
        """Write a FASTA header string from a ParsedHeader object."""
        return parsed_header.header


class UniprotParser:
    """Uniprot FASTA header parser."""

    header_pattern = re.compile(
        r"^(?P<db>\w+)\|(?P<id>[-\w]+)\|(?P<entry>\w+)\s+(?P<name>.*?)"
        r"(?:(\s+OS=(?P<OS>[^=]+))|"
        r"(\s+OX=(?P<OX>\d+))|"
        r"(\s+GN=(?P<GN>\S+))|"
        r"(\s+PE=(?P<PE>\d))|"
        r"(\s+SV=(?P<SV>\d+)))*\s*$"
    )

    field_names = {
        "db": "db",
        "id": "identifier",
        "entry": "entry_name",
        "name": "protein_name",
        "OS": "organism_name",
        "OX": "organism_identifier",
        "GN": "gene_name",
        "PE": "protein_existence",
        "SV": "sequence_version",
    }

    @classmethod
    def parse(cls, header: str) -> ParsedHeader:
        """Parse a FASTA header string into a ParsedHeader object."""
        match = cls.header_pattern.match(header)
        if match is None:
            raise ValueError(f"Header does not match the UniProt pattern: {header}")
        fields = match.groupdict()

        for key in ["OS", "OX", "GN", "PE", "SV"]:
            if fields[key] is None:
                del fields[key]
        fields = {cls.field_names[key]: value for key, value in fields.items()}

        return ParsedHeader(fields["identifier"], header, fields)

    @classmethod
    def write(cls, parsed_header: AbstractParsedHeader) -> str:
        """Write a FASTA header string from a ParsedHeader object."""
        fields = parsed_header.header_fields
        header_entries = [
            f"{fields['db']}|{fields['identifier']}|{fields['entry_name']}",
            f"{fields['protein_name']}",
        ]
        for key in ["OS", "OX", "GN", "PE", "SV"]:
            field_name = cls.field_names[key]
            if field_name not in fields:
                continue
            header_entries.append(f"{key}={fields[field_name]}")
        return " ".join(header_entries)


class UniprotLikeParser:
    """A tolerant FASTA header parser for UniProt like headers."""  
    field_pattern = re.compile(
        r"(?:(\s+OS=(?P<OS>[^=]+))|"
        r"(\s+OX=(?P<OX>\d+))|"
        r"(\s+GN=(?P<GN>\S+))|"
        r"(\s+PE=(?P<PE>\d))|"
        r"(\s+SV=(?P<SV>\d+)))*\s*$"        
    )

    tag_names = {
        "OS": "organism_name",
        "OX": "organism_identifier",
        "GN": "gene_name",
        "PE": "protein_existence",
        "SV": "sequence_version",
    }

    @classmethod
    def parse(cls, header: str) -> ParsedHeader:
        """Parse a FASTA header string into a ParsedHeader object."""
        split_header = header.split(maxsplit=1)
        try:
            db, _id, entry = split_header[0].split("|")
        except ValueError:
            raise ValueError(
                f"Header does not match the minimal UniProt like pattern: {header}"
            )
        fields = {"db": db, "identifier": _id, "entry_name": entry}

        if len(split_header) == 1:
            return ParsedHeader(fields["identifier"], header, fields)
        
        description = split_header[1]
        tag_positions = [description.find(f"{tag}=") for tag in cls.tag_names]
        matched_start = sorted([num for num in tag_positions if num >= 0])
        matched_end = matched_start[1:] + [len(description)]

        if not matched_start:  # Description contains only the protein name
            fields["protein_name"] = description
        elif matched_start[0] != 0:  # Description contains protein name and tag fields
            fields["protein_name"] = description[:matched_start[0]].rstrip()

        if matched_start:  # Header contains tag fields
            for start, end in zip(matched_start, matched_end):
                matched_field = description[start:end].rstrip().split("=", maxsplit=1)
                fields[matched_field[0]] = matched_field[1]
        
        for old_tag, new_tag in cls.tag_names.items():
            if old_tag in fields:
                fields[new_tag] = fields.pop(old_tag)

        return ParsedHeader(fields["identifier"], header, fields)

    @classmethod
    def write(cls, parsed_header: AbstractParsedHeader) -> str:
        """Write a FASTA header string from a ParsedHeader object."""
        fields = parsed_header.header_fields
        header_entries = [
            f"{fields['db']}|{fields['identifier']}|{fields['entry_name']}",
        ]
        if "protein_name" in fields:
            header_entries.append(f"{fields['protein_name']}")

        for key in ["OS", "OX", "GN", "PE", "SV"]:
            field_name = cls.field_names[key]
            if field_name not in fields:
                continue
            header_entries.append(f"{key}={fields[field_name]}")
        return " ".join(header_entries)


def register_parser(name: str, parser: HeaderParser):
    """Register a custom parser by name."""
    PARSER_REGISTRY[name] = parser


def get_parser(parser_name: str) -> HeaderParser:
    """Get a registered parser by name."""
    return PARSER_REGISTRY[parser_name]


PARSER_REGISTRY: dict[str, HeaderParser] = {
    "default": DefaultParser,
    "uniprot": UniprotParser,
    "uniprot_like": UniprotLikeParser,
}
