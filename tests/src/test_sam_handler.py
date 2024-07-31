from __future__ import annotations

from dataclasses import dataclass

import pytest

from csvtag.sam_handler import (
    _is_complete_overlapped,
    _is_non_microhomologic_overlapped,
    calculate_alignment_length,
    extract_alignment,
    extract_sqheaders,
    is_forward_strand,
    remove_overlapped_alignments,
)


def test_is_forward_strand():
    test_cases = [
        (0, True),  # Forward strand
        (16, False),  # Reverse strand
        (99, True),  # Forward strand
        (83, False),  # Reverse strand
        (147, False),  # Reverse strand
        (163, True),  # Forward strand
    ]

    for flag, expected in test_cases:
        assert is_forward_strand(flag) == expected, f"Failed for flag: {flag}"


def test_extract_sqheaders_simple():
    sam = [["@SQ", "SN:1", "LN:100"]]
    result = extract_sqheaders(sam)
    expected = {"1": 100}
    assert result == expected, f"Expected {expected}, but got {result}"


def test_extract_sqheaders_multiple():
    sam = [["@SQ", "SN:1", "LN:100", "SP:human"], ["@SQ", "SN:2", "LN:200", "SP:human"]]
    result = extract_sqheaders(sam)
    expected = {"1": 100, "2": 200}
    assert result == expected, f"Expected {expected}, but got {result}"


def test_extract_alignment():
    sam = [
        ["@SQ", "SN:1", "LN:100"],
        ["r001", "99", "chr1", "7", "255", "30M", "*", "0", "0", "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG", "*", "cs:Z::0"],
        ["r002", "0", "*", "0", "0", "*", "*", "0", "0", "*", "*", "cs:Z:1"],
    ]

    expected = [
        {
            "QNAME": "r001",
            "FLAG": 99,
            "RNAME": "chr1",
            "POS": 7,
            "MAPQ": 255,
            "CIGAR": "30M",
            "SEQ": "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG",
            "QUAL": "*",
            "CSTAG": ":0",
        }
    ]

    result = list(extract_alignment(sam))
    assert result == expected, f"Expected {expected}, but got {result}"


@pytest.mark.parametrize(
    "cigar, expected",
    [
        ("10M", 10),
        ("5M5I5M", 10),
        ("5M3D5M", 13),
        ("5M3N5M", 13),
        ("5M5S5M", 10),
    ],
)
def testcalculate_alignment_length(cigar, expected):
    result = calculate_alignment_length(cigar)
    assert result == expected, f"Expected {expected}, but got {result}"


#########################################
# _is_overlapped
#########################################


@dataclass
class OverlappedAlignment:
    prev_cigar: str
    curr_cigar: str
    prev_cstag: str = ""
    curr_cstag: str = ""
    prev_pos: int = 0
    curr_pos: int = 0


@pytest.mark.parametrize(
    "alignments_overlapped, expected",
    [
        (OverlappedAlignment("10M", "5M", prev_pos=1, curr_pos=1), True),
        (OverlappedAlignment("10M", "5M", prev_pos=1, curr_pos=6), True),
        (OverlappedAlignment("10M", "5M", prev_pos=1, curr_pos=7), False),
        (OverlappedAlignment("5M", "10M", prev_pos=1, curr_pos=1), False),
        (OverlappedAlignment("10M", "10M", prev_pos=1, curr_pos=1), True),
        (OverlappedAlignment("10M5I10M", "5M", prev_pos=1, curr_pos=11), True),
        (OverlappedAlignment("10M5I10M", "15M", prev_pos=1, curr_pos=11), False),
    ],
)
def test_is_complete_overlapped(alignments_overlapped, expected):
    result = _is_complete_overlapped(alignments_overlapped)
    assert result == expected, f"Expected {expected}, but got {result}"


@pytest.mark.parametrize(
    "alignments_overlapped, expected",
    [
        pytest.param(
            OverlappedAlignment(
                prev_cstag="=ACTG", curr_cstag="=ACT", prev_cigar="4M", curr_cigar="3M", prev_pos=1, curr_pos=1
            ),
            True,
            id="test_case_1",
        ),
        pytest.param(
            OverlappedAlignment(
                prev_cstag="=ACTG", curr_cstag="=ACT", prev_cigar="4M", curr_cigar="3M", prev_pos=1, curr_pos=6
            ),
            False,
            id="test_case_2",
        ),
        pytest.param(
            OverlappedAlignment(
                prev_cstag="=ACTG", curr_cstag="=ACT", prev_cigar="4M", curr_cigar="3M", prev_pos=1, curr_pos=7
            ),
            False,
            id="test_case_3",
        ),
        pytest.param(
            OverlappedAlignment(
                prev_cstag="=ACT", curr_cstag="=ACTG", prev_cigar="3M", curr_cigar="4M", prev_pos=1, curr_pos=1
            ),
            True,
            id="test_case_4",
        ),
        pytest.param(
            OverlappedAlignment(
                prev_cstag="=ACTG", curr_cstag="=ACTG", prev_cigar="4M", curr_cigar="4M", prev_pos=1, curr_pos=1
            ),
            True,
            id="test_case_5",
        ),
        pytest.param(
            OverlappedAlignment(
                prev_cstag="=ACTGACTGACTG",
                curr_cstag="=ACTGACTGACTG",
                prev_cigar="12M",
                curr_cigar="12M",
                prev_pos=1,
                curr_pos=11,
            ),
            True,
            id="test_case_6",
        ),
    ],
)
def test_is_non_microhomologic_overlapped(alignments_overlapped, expected):
    result = _is_non_microhomologic_overlapped(alignments_overlapped)
    assert result == expected, f"Expected {expected}, but got {result}"


@pytest.mark.parametrize(
    "alignments, expected",
    [
        pytest.param(
            [
                {"QNAME": "read1", "POS": 1, "CIGAR": "10M", "CSTAG": "=ACTGACTGAC"},
                {"QNAME": "read2", "POS": 1, "CIGAR": "10M", "CSTAG": "=TGCATGCATG"},
            ],
            [
                {"QNAME": "read1", "POS": 1, "CIGAR": "10M", "CSTAG": "=ACTGACTGAC"},
                {"QNAME": "read2", "POS": 1, "CIGAR": "10M", "CSTAG": "=TGCATGCATG"},
            ],
            id="no_duplicates",
        ),
        pytest.param(
            [
                {"QNAME": "read1", "POS": 1, "CIGAR": "5M", "CSTAG": "=ACTGA"},
                {"QNAME": "read1", "POS": 3, "CIGAR": "5M", "CSTAG": "=TGACT"},
            ],
            [{"QNAME": "read1", "POS": 1, "CIGAR": "5M", "CSTAG": "=ACTGA"}],
            id="overlapping_reads",
        ),
        pytest.param(
            [
                {"QNAME": "read1", "POS": 1, "CIGAR": "10M", "CSTAG": "=ACTGACTGAC"},
                {"QNAME": "read1", "POS": 1, "CIGAR": "5M", "CSTAG": "=ACTGA"},
            ],
            [{"QNAME": "read1", "POS": 1, "CIGAR": "10M", "CSTAG": "=ACTGACTGAC"}],
            id="contained_reads",
        ),
        pytest.param(
            [
                {"QNAME": "read1", "POS": 1, "CIGAR": "5M", "CSTAG": "=ACTGA"},
                {"QNAME": "read1", "POS": 3, "CIGAR": "5M", "CSTAG": "=TGACC"},
            ],
            [{"QNAME": "read1", "POS": 1, "CIGAR": "5M", "CSTAG": "=ACTGA"}],
            id="partially_overlapping_reads",
        ),
    ],
)
def test_remove_resequence(alignments, expected):
    result = list(remove_overlapped_alignments(iter(alignments)))
    assert result == expected, f"Expected {expected}, but got {result}"
