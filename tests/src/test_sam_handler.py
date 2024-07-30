from __future__ import annotations

import pytest

from csvtag.sam_handler import (
    _calculate_alignment_length,
    _is_complete_overlapped,
    extract_alignment,
    extract_sqheaders,
    remove_overlapped_alignments,
)


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
        ["r001", "99", "chr1", "7", "255", "30M", "*", "0", "0", "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG", "*", "cs:Z:0"],
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
def test_calculate_alignment_length(cigar, expected):
    result = _calculate_alignment_length(cigar)
    assert result == expected, f"Expected {expected}, but got {result}"


@pytest.mark.parametrize(
    "prev_cigar, prev_pos, curr_cigar, curr_pos, expected",
    [
        ("10M", 1, "5M", 1, True),
        ("10M", 1, "5M", 6, True),
        ("10M", 1, "5M", 7, False),
        ("5M", 1, "10M", 1, False),
        ("10M", 1, "10M", 1, True),
        ("10M5I10M", 1, "5M", 11, True),
        ("10M5I10M", 1, "15M", 11, False),
    ],
)
def test_is_complete_overlapped(prev_cigar, prev_pos, curr_cigar, curr_pos, expected):
    result = _is_complete_overlapped(prev_cigar, prev_pos, curr_cigar, curr_pos)
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
