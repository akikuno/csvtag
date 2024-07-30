from __future__ import annotations

from csvtag.sam_handler import extract_alignment, extract_sqheaders


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
