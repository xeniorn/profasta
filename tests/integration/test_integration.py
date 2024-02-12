import pytest
import profasta
import tempfile
from pathlib import Path


FASTA_PATH = Path(__file__).parent / "uniprot_hsapiens_10entries.fasta"


class CustomParser:
    @classmethod
    def parse(cls, header):
        return profasta.parser.ParsedHeader(header, header)

    @classmethod
    def write(self, parsed_header):
        return parsed_header.header


def test_built_in_uniprot_header_parser_read_write_header_roundtrip():
    UniprotParser = profasta.parser.get_parser("uniprot")

    with open(FASTA_PATH, "r") as file:
        for fasta_record in profasta.io.parse_fasta(file):
            parsed_header = UniprotParser.parse(fasta_record.header)
            written_header = UniprotParser.write(parsed_header)
            assert fasta_record.header == written_header


def test_custom_header_parser_read_write_header_roundtrip():
    profasta.parser.register_parser("custom_test_parser", CustomParser)
    PlainParser = profasta.parser.get_parser("custom_test_parser")

    with open(FASTA_PATH, "r") as file:
        for fasta_record in profasta.io.parse_fasta(file):
            parsed_header = PlainParser.parse(fasta_record.header)
            written_header = PlainParser.write(parsed_header)
            assert fasta_record.header == written_header


def test_protein_database_read_write_fasta_roundtrip(tmp_path):
    temp_fasta_path = tmp_path / "written_from_db.fasta"
    db = profasta.db.ProteinDatabase()
    db.add_fasta(FASTA_PATH, header_parser="default")
    db.write_fasta(temp_fasta_path, header_parser="default")

    db2 = profasta.db.ProteinDatabase()
    db2.add_fasta(temp_fasta_path, header_parser="default")

    for identifier in db:
        assert db[identifier].header == db2[identifier].header
        assert db[identifier].sequence == db2[identifier].sequence
