
import re

class Constants:
    class ProteinSequence:    
        terminator_symbol: str = '*'
        # protein sequence has at least 2 amino acids        
        valid_regex_string: str = '^[ACDEFGHIKLMNPQRSTVXY]{2,}$'
        # be gracious and allow typical alignment symbols and space stuff        
        ignore_symbols: list[str] = ['-', '+', '.', ' ', "\t"]        
            

def is_canonical_protein_sequence(input_sequence: str) -> bool:
    """Check whether the input sequence is a valid canonical representation.

    Canonical representation uses uppercase and does not include the terminator symbol.
    
    Args:
        input_sequence: The sequence to check.            

    Returns:
        True or False.
    """

    pattern: str = Constants.ProteinSequence.valid_regex_string
    result = re.fullmatch(pattern, input_sequence)

    is_valid = result is not None

    return is_valid

        
    

def make_canonical_sequence(input_sequence: str) -> str:
    """Take an arbitarily formatted string and return a canonical sequence representation.

    Canonical representation uses uppercase and does not include the terminator symbol.
    
    Args:
        input_sequence: The sequence to canonicalize.            

    Returns:
        Canonical sequence.

    Throws:
        Exception, when input cannot be unambiguously translated into a valid protein sequence
    
    """ 

    seq: str = input_sequence.upper()

    # remove whitespace if exists. Do this first!
    seq = ''.join(seq.split())    

    terminator: str = Constants.ProteinSequence.terminator_symbol
    
    # remove terminator from the end, even if there are multiple
    while seq.endswith(terminator):
        seq = seq.removesuffix(terminator)
    
    # really weird situation but allow it nevertheless if it's at the very beginning too
    while seq.startswith(terminator):
        seq = seq.removeprefix(terminator)

    # these are ok to be inside and are just cut out and fully ignored
    ignore_symbols = Constants.ProteinSequence.ignore_symbols
    for symbol in ignore_symbols:
        seq = seq.replace(symbol, '')   

    # e.g. has a terminator symbol in the middle, or otherwise other disallowed symbols
    if not is_canonical_protein_sequence(seq):
        message = f"Input cannot be parsed into a canonical sequence: {input_sequence}\nTried until: {seq}"
        raise Exception(message)
    
    return seq
