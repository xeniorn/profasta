import io
import pytest

import profasta.io


def test_parse_fasta_file():
    file_buffer = io.StringIO(f">H1\nMKKK\n>H2\nMAAA")
    fasta_parser = profasta.io.parse_fasta(file_buffer)
    record = next(fasta_parser)

    assert record.header == "H1"
    assert record.sequence == "MKKK"


def test_write_fasta_file():
    fasta_records = [
        profasta.io.FastaRecord(header="Header1", sequence="ACGT"),
        profasta.io.FastaRecord(header="Header2", sequence="TGCA"),
    ]

    file_buffer = io.StringIO()
    profasta.io.write_fasta(file_buffer, fasta_records)

    file_buffer.seek(0)
    file_contents = file_buffer.read()

    expected_output = ">Header1\nACGT\n>Header2\nTGCA\n"
    assert file_contents == expected_output


class TestMakeRecordString:
    @pytest.fixture
    def setup(self):
        self.header = "Header"
        self.sequence = "ACGTACGTACGTACGT"

    def test_with_no_line_breaks(self, setup):
        expected = ">Header\nACGTACGTACGTACGT"
        assert profasta.io.make_record_string(self.header, self.sequence, line_width=-1) == expected  # fmt: skip

    def test_with_line_width_5(self, setup):
        expected = ">Header\nACGTA\nCGTAC\nGTACG\nT"
        assert profasta.io.make_record_string(self.header, self.sequence, line_width=5) == expected  # fmt: skip

    def test_with_line_width_longer_than_the_sequence(self, setup):
        expected = ">Header\nACGTACGTACGTACGT"
        assert profasta.io.make_record_string(self.header, self.sequence, line_width=99) == expected  # fmt: skip
