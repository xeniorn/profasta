import pytest

import profasta.parser


class TestUniprotLikeParser:
    def test_parse_full_uniprot_header(self):
        header = "sp|O75385|ULK1_HUMAN Serine/threonine-protein kinase ULK1 OS=Homo sapiens OX=9606 GN=ULK1 PE=1 SV=2"
        expected_fields = {
            "db": "sp",
            "identifier": "O75385",
            "entry_name": "ULK1_HUMAN",
            "protein_name": "Serine/threonine-protein kinase ULK1",
            "organism_name": "Homo sapiens",
            "organism_identifier": "9606",
            "gene_name": "ULK1",
            "protein_existence": "1",
            "sequence_version": "2",
        }
        parsed_header = profasta.parser.UniprotLikeParser.parse(header)
        assert parsed_header.header_fields == expected_fields

    def test_parse_partial_uniprot_header_with_some_tag_fields_missing(self):
        header = (
            "sp|O75385|ULK1_HUMAN Serine/threonine-protein kinase ULK1 OX=9606 GN=ULK1"
        )
        expected_fields = {
            "db": "sp",
            "identifier": "O75385",
            "entry_name": "ULK1_HUMAN",
            "protein_name": "Serine/threonine-protein kinase ULK1",
            "organism_identifier": "9606",
            "gene_name": "ULK1",
        }
        parsed_header = profasta.parser.UniprotLikeParser.parse(header)
        assert parsed_header.header_fields == expected_fields

    def test_parse_partial_uniprot_header_with_protein_name_missing(self):
        header = "sp|O75385|ULK1_HUMAN OS=Homo sapiens OX=9606 GN=ULK1 PE=1 SV=2"
        expected_fields = {
            "db": "sp",
            "identifier": "O75385",
            "entry_name": "ULK1_HUMAN",
            "organism_name": "Homo sapiens",
            "organism_identifier": "9606",
            "gene_name": "ULK1",
            "protein_existence": "1",
            "sequence_version": "2",
        }
        parsed_header = profasta.parser.UniprotLikeParser.parse(header)
        assert parsed_header.header_fields == expected_fields

    def test_parse_partial_uniprot_header_with_no_tag_fields(self):
        header = "sp|O75385|ULK1_HUMAN Serine/threonine-protein kinase ULK1"
        expected_fields = {
            "db": "sp",
            "identifier": "O75385",
            "entry_name": "ULK1_HUMAN",
            "protein_name": "Serine/threonine-protein kinase ULK1",
        }

        parsed_header = profasta.parser.UniprotLikeParser.parse(header)
        assert parsed_header.header_fields == expected_fields

    def test_parse_minimal_uniprot_header_with_no_description(self):
        header = "sp|O75385|ULK1_HUMAN"
        expected_fields = {
            "db": "sp",
            "identifier": "O75385",
            "entry_name": "ULK1_HUMAN",
        }

        parsed_header = profasta.parser.UniprotLikeParser.parse(header)
        assert parsed_header.header_fields == expected_fields


def test_decoy_writer():
    header = "sp|O75385|ULK1_HUMAN"
    parsed_header = profasta.parser.ParsedHeader("", header, {})
    decoy_header = profasta.parser.DecoyWriter.write(parsed_header)
    assert decoy_header == "rev_" + header
