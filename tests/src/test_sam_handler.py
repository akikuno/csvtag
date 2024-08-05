from __future__ import annotations

import pytest

from csvtag.sam_handler import (
    calculate_alignment_length,
    extract_alignment,
    extract_sqheaders,
    is_forward_strand,
    trim_softclip,
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


@pytest.mark.parametrize(
    "qual, cigar, expected",
    [
        ("ABCDEFGHI", "2S7M", "CDEFGHI"),  # Remove softclip at the beginning
        ("ABCDEFGHI", "7M2S", "ABCDEFG"),  # Remove softclip at the end
        ("ABCDEFGHI", "2S5M2S", "CDEFG"),  # Remove softclip at both ends
        ("ABCDEFGHI", "9M", "ABCDEFGHI"),  # No softclip
        ("ABCDEFGHI", "1S8M", "BCDEFGHI"),  # Remove 1S at the beginning
        ("ABCDEFGHI", "8M1S", "ABCDEFGH"),  # Remove 1S at the end
    ],
)
def test_trim_softclip(qual, cigar, expected):
    assert trim_softclip(qual, cigar) == expected
