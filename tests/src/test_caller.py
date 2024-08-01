from __future__ import annotations

from pathlib import Path

import pytest

from csvtag.caller import _is_second_strand_different, _is_within_bases, _padding_n, call_csvtag


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


@pytest.mark.parametrize(
    "cs_tag, length, side, expected",
    [
        ("=ACGT", 2, "left", "=NNACGT"),
        ("*AC=ACGT", 2, "left", "=NN*AC=ACGT"),
        ("=ACGT", 2, "right", "=ACGTNN"),
        ("=A*TC=GT", 2, "right", "=A*TC=GTNN"),
        ("=A*TC", 2, "right", "=A*TC=NN"),
    ],
)
def test_padding_n(cs_tag, length, side, expected):
    result = _padding_n(cs_tag, length, side)
    assert result == expected, f"Expected {expected}, but got {result}"


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


# TODO: Add inversion tests


def test_call_csvtag_three_alignments_with_inv():
    path_sam = Path("tests/data/three_alignments_witn_inv.sam")
    result = list(call_csvtag(path_sam))
    expected = [
        {"QNAME": "read1", "CSVTAG": "=AAAAANNNNN=aa*ag=aa=NNNNNGGGGG", "POS": 1, "RNAME": "ref"},
    ]
    assert result == expected, f"Expected {expected}, but got {result}"


def test_call_csvtag_four_alignments():
    path_sam = Path("tests/data/four_alignments.sam")
    result = list(call_csvtag(path_sam))
    result.sort(key=lambda x: [x["QNAME"], x["POS"]])
    expected = [
        {"QNAME": "read1", "CSVTAG": "=AAAAANNNNN=aa*ag=aa=NNNNNGGGGG", "POS": 1, "RNAME": "ref"},
        {"QNAME": "read1", "CSVTAG": "=ACGT", "POS": 201, "RNAME": "ref"},
        {"QNAME": "read2", "CSVTAG": "=ACGT", "POS": 1, "RNAME": "ref"},
        {"QNAME": "read2", "CSVTAG": "=AAAAANNNNN=aa*ag=aa=NNNNNGGGGG", "POS": 201, "RNAME": "ref"},
    ]
    expected.sort(key=lambda x: [x["QNAME"], x["POS"]])
    assert result == expected, f"Expected {expected}, but got {result}"
