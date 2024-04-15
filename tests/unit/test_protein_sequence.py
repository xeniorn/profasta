import pytest

import profasta.protein_sequence


class TestProteinSequenceParsers:
    
    @pytest.mark.parametrize(
        "sequence, expected_result",
        [
            ("ACDEFGHIKLMNPQRSTVXY", "ACDEFGHIKLMNPQRSTVXY"),
            ("ACD ACD", "ACDACD"),
            ("Xx", "XX"),
            ("acdefghiklmnpqrstvxy", "ACDEFGHIKLMNPQRSTVXY"),            
            ("acdefghiklmnpqrsTVXY", "ACDEFGHIKLMNPQRSTVXY"),            
            ("ACDEfghiklmnpqrstvxy", "ACDEFGHIKLMNPQRSTVXY"),                                    
            ("acDEfghiKLmnPQRSTvXy", "ACDEFGHIKLMNPQRSTVXY"),                                                
            ("xxATx", "XXATX")            
        ],
    )    
    def test_canonicalize_success_valid_sequences(self, sequence, expected_result):                
        result: str = profasta.protein_sequence.make_canonical_sequence(sequence)
        assert result == expected_result


    @pytest.mark.parametrize(
        "sequence, expected_result",
        [
            ("AG.GA", "AGGA"),
            ("AGR-GA", "AGRGA"),
            ("AHGGA*", "AHGGA"),
            ("*AQGGA", "AQGGA"),
            ("acdefghik lmnpqrstvxy", "ACDEFGHIKLMNPQRSTVXY"),                                    
            ("acdefghik" + "\t" + "lmnpqrstvxy", "ACDEFGHIKLMNPQRSTVXY"),                                    
            ("**aC-d.e++fgHIk" + "\t \t" + "l..m---n" + "\t" + "pq rs tv XY* **", "ACDEFGHIKLMNPQRSTVXY"),                                                
        ],
    )    
    def test_canonicalize_success_odd_but_unambiguous_sequences(self, sequence, expected_result):                
        result: str = profasta.protein_sequence.make_canonical_sequence(sequence)
        assert result == expected_result        

    @pytest.mark.parametrize(
        "sequence",
        [
            "",
            "A",
            "AGC1",
            "1QVT",
            "TCRAQ75VNMAGGC",                        
            "AG#GA",
            "AG?GA",            
            "AG*GA",
            "AZ",
            "IOT",
            "JURAJ",
            "BAHEL"                                                                                     
        ],
    )    
    def test_canonicalize_fail_invalid_sequences(self, sequence):
        with pytest.raises(Exception):
            result: str = profasta.protein_sequence.make_canonical_sequence(sequence)
        
        
def tryme():
    sequence = "**aC-d.e++fgHIk" + "\t \t" + "l..m---n" + "\t" + "pq rs tv XY* **"
    profasta.protein_sequence.make_canonical_sequence(sequence)