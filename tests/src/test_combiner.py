from __future__ import annotations

import pytest

from csvtag.combiner import _padding_n, combine_neighboring_csv_tags, combine_splitted_tags


@pytest.mark.parametrize(
    "splitted_csv_tag, expected",
    [
        (iter(["=A", "+T|+T|+T|=C", "=C", "-A", "-A", "=T", "*AG", "=T", "=T"]), "=A+TTT=CC-AA=T*AG=TT"),
        (iter(["=A", "=A", "=A", "=N", "=N", "=N", "=N", "=N", "=C", "=C", "=A"]), "=AAANNNNNCCA"),
        (iter(["*AG", "*T", "*C", "*G"]), "*AG*T*C*G"),
        (iter(["=G", "-A", "+T", "=C"]), "=G-A+T=C"),
    ],
)
def test_combine_splitted_tags(splitted_csv_tag, expected):
    result = combine_splitted_tags(splitted_csv_tag)
    assert result == expected, f"Expected {expected}, but got {result}"


###########################################################
# combine_neighboring_csv_tags
###########################################################


@pytest.mark.parametrize(
    "csv_tag, n_length, side, expected",
    [
        ("=ACGT", 2, "left", "=NNACGT"),
        ("*AC=ACGT", 2, "left", "=NN*AC=ACGT"),
        ("=ACGT", 2, "right", "=ACGTNN"),
        ("=A*TC=GT", 2, "right", "=A*TC=GTNN"),
        ("=A*TC", 2, "right", "=A*TC=NN"),
    ],
)
def test_padding_n(csv_tag, n_length, side, expected):
    result = _padding_n(csv_tag, n_length, side)
    assert result == expected, f"Expected {expected}, but got {result}"


@pytest.mark.parametrize(
    "csv_tags, positions, distance, expected",
    [
        (["=AA", "=tt", "=CC"], [1, 5, 10], 50, ["=AANN=tt=NNNCC"]),
        # (["=AA", "=tt", "=CC"], [1, 5, 100], 50, ["=AANN=tt", "=CC"]),
        # (
        #     [
        #         "=GGGGTATGTGGCTGCGTGGTCAAATGTGCGGCATACGTATTTGCTCGGCGTGCTTGCTCTCACGAACTTGACCTGGAGATCAAGGAGATGTTTCTTGTCCAACTGGACAGCGCTTCAACGGAATGGATCTACGTTACAGCCTGCATAAAGAAAACGGAGTTGCCGAGGACGAAAGCGACTTTAGGTTCTGTCCGTTGTCTTTGGCGGAAAACTTCCACTCAGGAAGCAGACACTGATTGACACGGTTTAGCA",
        #         "=gcacttgatcgttttgctgtagaaaaaacttaataaacagaatgccgatgaaggcactactgtactaatagggccgggctacatgttaactacaaggctataacctattgatgacccggtccatacataacttggtatcgtgcatgtagcgttcaagggctatagcaattccgacggaaatccattggggtaacgccttagaataatatactggcctatcgcaacacaaccacctctgccgtgtaatccgagggtggccacgacaatcgaaggtatggtcgaccgttgtaggtaattctaggcgatgaggggtccttctttcataaattttcttcaggacattgttcacgtaaactaccaggattaaccgtcgtagtgagcccgcttggttttgggaactcgtgtcttaaattcgtctccgattagcgcactatactatactttaatcccagacataccgatattaaaccactcaatttgacctaatcctcaaaccttctgc",
        #         "=TGCACTTCCACAGAGCGCGGTAGAGACTCATCCACCCGGCAGCTCTGTAATAGGGACTAAAAAAGTGATGATAATCATGAGTGCCGCGTTATGGTGGTGTCGGATCAGAGCGGTCTTACGACCAGTCGTATGCCTTCTCGAGTTCCGTCCGGTTAAGCGTGACAGTCCCAGTGAACCCACAAACCGTGATGGCTGTCCTTGGAGTCATACGCAAGAAGGATGGTCTCCAGACACCGGCGCACCAGTTTTCAC",
        #     ],
        #     [1, 250, 749],
        #     50,
        #     [
        #         "=GGGGTATGTGGCTGCGTGGTCAAATGTGCGGCATACGTATTTGCTCGGCGTGCTTGCTCTCACGAACTTGACCTGGAGATCAAGGAGATGTTTCTTGTCCAACTGGACAGCGCTTCAACGGAATGGATCTACGTTACAGCCTGCATAAAGAAAACGGAGTTGCCGAGGACGAAAGCGACTTTAGGTTCTGTCCGTTGTCTTTGGCGGAAAACTTCCACTCAGGAAGCAGACACTGATTGACACGGTTTA=gcacttgatcgttttgctgtagaaaaaacttaataaacagaatgccgatgaaggcactactgtactaatagggccgggctacatgttaactacaaggctataacctattgatgacccggtccatacataacttggtatcgtgcatgtagcgttcaagggctatagcaattccgacggaaatccattggggtaacgccttagaataatatactggcctatcgcaacacaaccacctctgccgtgtaatccgagggtggccacgacaatcgaaggtatggtcgaccgttgtaggtaattctaggcgatgaggggtccttctttcataaattttcttcaggacattgttcacgtaaactaccaggattaaccgtcgtagtgagcccgcttggttttgggaactcgtgtcttaaattcgtctccgattagcgcactatactatactttaatcccagacataccgatattaaaccactcaatttgacctaatcctcaaaccttc=TGCACTTCCACAGAGCGCGGTAGAGACTCATCCACCCGGCAGCTCTGTAATAGGGACTAAAAAAGTGATGATAATCATGAGTGCCGCGTTATGGTGGTGTCGGATCAGAGCGGTCTTACGACCAGTCGTATGCCTTCTCGAGTTCCGTCCGGTTAAGCGTGACAGTCCCAGTGAACCCACAAACCGTGATGGCTGTCCTTGGAGTCATACGCAAGAAGGATGGTCTCCAGACACCGGCGCACCAGTTTTCAC"
        #     ],
        # ),
    ],
)
def test_combine_neighboring_csv_tags(csv_tags, positions, distance, expected):
    result = list(combine_neighboring_csv_tags(csv_tags, positions, distance))
    assert result == expected, f"Expected {expected}, but got {result}"
