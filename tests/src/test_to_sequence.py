from __future__ import annotations

import pytest
from csvtag.to_sequence import to_sequence


@pytest.mark.parametrize(
    "csv_tag, expected",
    [
        ("=AA=aa*ga=a=AA", "AAaaaaAA"),  # Basic CSV tag
        ("=AA=aa+gg=aa=AA", "AAaaggaaAA"),  # With insertion
        ("=AA=aa-gg=aa=AA", "AAaaaaAA"),  # With deletion
        ("=AAGG*CT=AT", "AAGGTAT"),  # All uppercase CSV tag
        ("=tt", "tt"),  # All lowercase CSV tag
        ("=A=a=G", "AaG"),  # Mixed uppercase and lowercase CSV tag
        ("", ""),  # Empty CSV tag
        (
            "=GGGGTATGTGGCTGCGTGGTCAAATGTGCGGCATACGTATTTGCTCGGCGTGCTTGCTCTCACGAACTTGACCTGGAGATCAAGGAGATGTTTCTTGTCCAACTGGACAGCGCTTCAACGGAATGGATCTACGTTACAGCCTGCATAAAGAAAACGGAGTTGCCGAGGACGAAAGCGACTTTAGGTTCTGTCCGTTGTCTTTGGCGGAAAACTTCCACTCAGGAAGCAGACACTGATTGACACGGTTTAGCA=cttgatcgttttgctgtagaaaaaacttaataaacagaatgccgatgaaggcactactgtactaatagggccgggctacatgttaactacaaggctataacctattgatgacccggtccatacataacttggtatcgtgcatgtagcgttcaagggctatagcaattccgacggaaatccattggggtaacgccttagaataatatactggcctatcgcaacacaaccacctctgccgtgtaatccgagggtggccacgacaatcgaaggtatggtcgaccgttgtaggtaattctaggcgatgaggggtccttctttcataaattttcttcaggacattgttcacgtaaactaccaggattaaccgtcgtagtgagcccgcttggttttgggaactcgtgtcttaaattcgtctccgattagcgcactatactatactttaatcccagacataccgatattaaaccactcaatttgacctaatcctcaaaccttctgc=ACTTCCACAGAGCGCGGTAGAGACTCATCCACCCGGCAGCTCTGTAATAGGGACTAAAAAAGTGATGATAATCATGAGTGCCGCGTTATGGTGGTGTCGGATCAGAGCGGTCTTACGACCAGTCGTATGCCTTCTCGAGTTCCGTCCGGTTAAGCGTGACAGTCCCAGTGAACCCACAAACCGTGATGGCTGTCCTTGGAGTCATACGCAAGAAGGATGGTCTCCAGACACCGGCGCACCAGTTTTCAC",
            "GGGGTATGTGGCTGCGTGGTCAAATGTGCGGCATACGTATTTGCTCGGCGTGCTTGCTCTCACGAACTTGACCTGGAGATCAAGGAGATGTTTCTTGTCCAACTGGACAGCGCTTCAACGGAATGGATCTACGTTACAGCCTGCATAAAGAAAACGGAGTTGCCGAGGACGAAAGCGACTTTAGGTTCTGTCCGTTGTCTTTGGCGGAAAACTTCCACTCAGGAAGCAGACACTGATTGACACGGTTTAGCActtgatcgttttgctgtagaaaaaacttaataaacagaatgccgatgaaggcactactgtactaatagggccgggctacatgttaactacaaggctataacctattgatgacccggtccatacataacttggtatcgtgcatgtagcgttcaagggctatagcaattccgacggaaatccattggggtaacgccttagaataatatactggcctatcgcaacacaaccacctctgccgtgtaatccgagggtggccacgacaatcgaaggtatggtcgaccgttgtaggtaattctaggcgatgaggggtccttctttcataaattttcttcaggacattgttcacgtaaactaccaggattaaccgtcgtagtgagcccgcttggttttgggaactcgtgtcttaaattcgtctccgattagcgcactatactatactttaatcccagacataccgatattaaaccactcaatttgacctaatcctcaaaccttctgcACTTCCACAGAGCGCGGTAGAGACTCATCCACCCGGCAGCTCTGTAATAGGGACTAAAAAAGTGATGATAATCATGAGTGCCGCGTTATGGTGGTGTCGGATCAGAGCGGTCTTACGACCAGTCGTATGCCTTCTCGAGTTCCGTCCGGTTAAGCGTGACAGTCCCAGTGAACCCACAAACCGTGATGGCTGTCCTTGGAGTCATACGCAAGAAGGATGGTCTCCAGACACCGGCGCACCAGTTTTCAC",
        ),
    ],
)
def test_to_sequence(csv_tag: str, expected: str):
    result = to_sequence(csv_tag)
    assert result == expected, f"Expected {expected}, but got {result}"
