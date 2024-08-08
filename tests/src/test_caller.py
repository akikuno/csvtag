from __future__ import annotations

from pathlib import Path

import pytest
from csvtag.caller import _is_second_strand_different, _is_within_bases, call_csvtag


@pytest.mark.parametrize(
    "first_flag, second_flag, third_flag, expected",
    [
        (0, 16, 0, True),  # first and third are forward, second is reverse
        (16, 0, 16, True),  # first and third are reverse, second is forward
        (0, 0, 0, False),  # all forward
        (16, 16, 16, False),  # all reverse
        (0, 0, 16, False),  # first and second are forward, third is reverse
        (16, 0, 0, False),  # first is reverse, second and third are forward
    ],
)
def test_is_second_strand_different(first_flag, second_flag, third_flag, expected):
    assert _is_second_strand_different(first_flag, second_flag, third_flag) == expected


@pytest.mark.parametrize(
    "first_end, second_pos, second_end, third_pos, expected",
    [
        (5, 10, 15, 20, True),
        (10, 100, 125, 135, False),  # second_pos - first_end = 11 => False
        (0, 5, 10, 100, False),
    ],
)
def test_is_within_bases(first_end, second_pos, second_end, third_pos, expected):
    assert _is_within_bases(first_end, second_pos, second_end, third_pos, base_num=10) == expected


def test_call_csvtag_one_alignment():
    path_sam = Path("tests/data/one_alignment.sam")
    result = list(call_csvtag(path_sam))
    expected = [{"QNAME": "read1", "CSVTAG": "=AAAAA", "POS": 1, "RNAME": "ref"}]
    assert result == expected, f"Expected {expected}, but got {result}"


def test_call_csvtag_two_alignments():
    path_sam = Path("tests/data/two_alignments.sam")
    result = list(call_csvtag(path_sam))
    expected = [
        {"QNAME": "read1", "CSVTAG": "=AAAAA", "POS": 1, "RNAME": "ref"},
        {"QNAME": "read1", "CSVTAG": "=AA*AG=AA", "POS": 100, "RNAME": "ref"},
    ]
    assert result == expected, f"Expected {expected}, but got {result}"


def test_call_csvtag_three_alignments_with_inv():
    path_sam = Path("tests/data/three_alignments_witn_inv.sam")
    result = list(call_csvtag(path_sam))
    expected = [
        {"QNAME": "read1", "RNAME": "ref", "POS": 1, "CSVTAG": "=AAAAA"},
        {"QNAME": "read1", "RNAME": "ref", "POS": 11, "CSVTAG": "=aa*ag=aa"},
        {"QNAME": "read1", "RNAME": "ref", "POS": 21, "CSVTAG": "=GGGGG"},
    ]
    assert result == expected, f"Expected {expected}, but got {result}"


def test_call_csvtag_four_alignments():
    path_sam = Path("tests/data/four_alignments.sam")
    result = list(call_csvtag(path_sam))
    result.sort(key=lambda x: [x["QNAME"], x["RNAME"], x["POS"]])
    expected = [
        {"QNAME": "read1", "RNAME": "ref", "POS": 1, "CSVTAG": "=AAAAA"},
        {"QNAME": "read1", "RNAME": "ref", "POS": 11, "CSVTAG": "=aa*ag=aa"},
        {"QNAME": "read1", "RNAME": "ref", "POS": 21, "CSVTAG": "=GGGGG"},
        {"QNAME": "read1", "RNAME": "ref", "POS": 201, "CSVTAG": "=ACGT"},
        {"QNAME": "read2", "RNAME": "ref", "POS": 1, "CSVTAG": "=ACGT"},
        {"QNAME": "read2", "RNAME": "ref", "POS": 201, "CSVTAG": "=AAAAA"},
        {"QNAME": "read2", "RNAME": "ref", "POS": 211, "CSVTAG": "=aa*ag=aa"},
        {"QNAME": "read2", "RNAME": "ref", "POS": 221, "CSVTAG": "=GGGGG"},
    ]
    expected.sort(key=lambda x: [x["QNAME"], x["POS"]])
    assert result == expected, f"Expected {expected}, but got {result}"


def test_call_csvtag_inversion_sr_simulated():
    path_sam = Path("tests/data/inversion_sr_simulated.sam")
    result = list(call_csvtag(path_sam))
    expected = [
        {
            "CSVTAG": "=GGGGTATGTGGCTGCGTGGTCAAATGTGCGGCATACGTATTTGCTCGGCGTGCTTGCTCTCACGAACTTGACCTGGAGATCAAGGAGATGTTTCTTGTCCAACTGGACAGCGCTTCAACGGAATGGATCTACGTTACAGCCTGCATAAAGAAAACGGAGTTGCCGAGGACGAAAGCGACTTTAGGTTCTGTCCGTTGTCTTTGGCGGAAAACTTCCACTCAGGAAGCAGACACTGATTGACACGGTTTAGCA",
            "POS": 1,
            "QNAME": "read1",
            "RNAME": "control",
        },
        {
            "CSVTAG": "=gcacttgatcgttttgctgtagaaaaaacttaataaacagaatgccgatgaaggcactactgtactaatagggccgggctacatgttaactacaaggctataacctattgatgacccggtccatacataacttggtatcgtgcatgtagcgttcaagggctatagcaattccgacggaaatccattggggtaacgccttagaataatatactggcctatcgcaacacaaccacctctgccgtgtaatccgagggtggccacgacaatcgaaggtatggtcgaccgttgtaggtaattctaggcgatgaggggtccttctttcataaattttcttcaggacattgttcacgtaaactaccaggattaaccgtcgtagtgagcccgcttggttttgggaactcgtgtcttaaattcgtctccgattagcgcactatactatactttaatcccagacataccgatattaaaccactcaatttgacctaatcctcaaaccttctgc",
            "POS": 250,
            "QNAME": "read1",
            "RNAME": "control",
        },
        {
            "CSVTAG": "=TGCACTTCCACAGAGCGCGGTAGAGACTCATCCACCCGGCAGCTCTGTAATAGGGACTAAAAAAGTGATGATAATCATGAGTGCCGCGTTATGGTGGTGTCGGATCAGAGCGGTCTTACGACCAGTCGTATGCCTTCTCGAGTTCCGTCCGGTTAAGCGTGACAGTCCCAGTGAACCCACAAACCGTGATGGCTGTCCTTGGAGTCATACGCAAGAAGGATGGTCTCCAGACACCGGCGCACCAGTTTTCAC",
            "POS": 749,
            "QNAME": "read1",
            "RNAME": "control",
        },
    ]
    result.sort(key=lambda x: [x["QNAME"], x["RNAME"], x["POS"]])
    expected.sort(key=lambda x: [x["QNAME"], x["RNAME"], x["POS"]])
    assert result == expected, f"Expected {expected}, but got {result}"
