from __future__ import annotations

import pytest

from csvtag.microhomology_handler import _get_length_of_microhomology, _trim_microhomology, _trim_softclip


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
    assert _trim_softclip(qual, cigar) == expected


@pytest.mark.parametrize(
    "curr_sequence, next_sequence, curr_qual, next_qual, expected",
    [
        ("AAAA", "GGGG", "!!!!", "!!!!", 0),  # No overlap
        ("AAGG", "GGCC", "!!!!", "!!!!", 2),  # 2 bases overlap
        ("AGCT", "AGCT", "!!!!", "!!!!", 4),  # "AGCT" overlap
        ("TAGCT", "AGCT", "!!!!!", "!!!!", 4),  # "AGCT" overlap
        ("TAGCT", "AGCTA", "!!!!!", "!!!!!", 4),  # "AGCT" overlap
        ("AAAA", "AAAA", "!!!!", "@@@@", 0),  # No overlap due to mismatch at quality
        ("AAAA", "AAAA", "!!!!", "!!@@", 2),  # "AA" overlap
    ],
)
def test_get_length_of_microhomology(curr_sequence, next_sequence, curr_qual, next_qual, expected):
    assert _get_length_of_microhomology(curr_sequence, next_sequence, curr_qual, next_qual) == expected


@pytest.mark.parametrize(
    "alignments, expected",
    [
        (
            [
                {"QNAME": "read1", "CSTAG": "=AAATTT", "CIGAR": "6M", "QUAL": "!!!!!!"},
                {"QNAME": "read1", "CSTAG": "=TTTCCC", "CIGAR": "6M", "QUAL": "!!!!!!"},
            ],
            [
                {"QNAME": "read1", "CSTAG": "=AAATTT", "CIGAR": "6M", "QUAL": "!!!!!!"},
                {"QNAME": "read1", "CSTAG": "=CCC", "CIGAR": "6M", "QUAL": "!!!!!!"},
            ],
        ),
        (
            [
                {"QNAME": "read1", "CSTAG": "=AAATTTGGG", "CIGAR": "9M", "QUAL": "!!!!!!!!!"},
                {"QNAME": "read1", "CSTAG": "=GGGCCC", "CIGAR": "6M", "QUAL": "!!!!!!"},
            ],
            [
                {"QNAME": "read1", "CSTAG": "=AAATTTGGG", "CIGAR": "9M", "QUAL": "!!!!!!!!!"},
                {"QNAME": "read1", "CSTAG": "=CCC", "CIGAR": "6M", "QUAL": "!!!!!!"},
            ],
        ),
        (
            [
                {"QNAME": "read1", "CSTAG": "=AAATTT", "CIGAR": "6M", "QUAL": "!!!!!!"},
                {"QNAME": "read1", "CSTAG": "=CCCGGG", "CIGAR": "6M", "QUAL": "!!!!!!"},
            ],
            [
                {"QNAME": "read1", "CSTAG": "=AAATTT", "CIGAR": "6M", "QUAL": "!!!!!!"},
                {"QNAME": "read1", "CSTAG": "=CCCGGG", "CIGAR": "6M", "QUAL": "!!!!!!"},
            ],
        ),
    ],
)
def test_trim_microhomology(alignments, expected):
    result = _trim_microhomology(alignments)
    assert result == expected, f"Expected {expected}, but got {result}"
