from __future__ import annotations

import pytest

from csvtag.microhomology_handler import _get_length_of_microhomology, _trim_microhomology, trim_microhomology

# @pytest.mark.parametrize(
#     "curr_sequence, next_sequence, curr_qual, next_qual, expected",
#     [
#         ("AAAA", "GGGG", "!!!!", "!!!!", 0),  # No overlap
#         ("AAGG", "GGCC", "!!!!", "!!!!", 2),  # 2 bases overlap
#         ("AGCT", "AGCT", "!!!!", "!!!!", 4),  # "AGCT" overlap
#         ("TAGCT", "AGCT", "!!!!!", "!!!!", 4),  # "AGCT" overlap
#         ("TAGCT", "AGCTA", "!!!!!", "!!!!!", 4),  # "AGCT" overlap
#         ("AAAA", "AAAA", "!!!!", "@@@@", 0),  # No overlap due to mismatch at quality
#         ("AAAA", "AAAA", "!!!!", "!!@@", 2),  # "AA" overlap
#     ],
# )
# def test_get_length_of_microhomology(curr_sequence, next_sequence, curr_qual, next_qual, expected):
#     assert _get_length_of_microhomology(curr_sequence, next_sequence, curr_qual, next_qual) == expected


@pytest.mark.parametrize(
    "curr_sequence, next_sequence, expected",
    [
        ("AAAA", "GGGG", 0),  # No overlap
        ("AAGG", "GGCC", 2),  # 2 bases overlap
        ("AGCT", "AGCT", 4),  # "AGCT" overlap
        ("TAGCT", "AGCT", 4),  # "AGCT" overlap
        ("TAGCT", "AGCTA", 4),  # "AGCT" overlap
        ("AAAA", "AAAA", 4),  # No overlap due to mismatch at quality
    ],
)
def test_get_length_of_microhomology(curr_sequence, next_sequence, expected):
    assert _get_length_of_microhomology(curr_sequence, next_sequence) == expected


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
def test__trim_microhomology(alignments, expected):
    result = _trim_microhomology(alignments)
    assert result == expected, f"Expected {expected}, but got {result}"


###########################################################
# Handle csv tag
###########################################################


@pytest.mark.parametrize(
    "csv_tags, expected",
    [
        (["=AAATTT", "=TTTCCC"], ["=AAATTT", "=CCC"]),
        (["=AAATTT", "=TTTCCC", "=CCGG"], ["=AAATTT", "=CCC", "=GG"]),
        (["=AAATTT", "=tttccc"], ["=AAATTT", "=ccc"]),
        (["=aaattt", "=TTTCCC"], ["=aaattt", "=CCC"]),
        (["=A", "=C", "=G", "=T"], ["=A", "=C", "=G", "=T"]),
    ],
)
def test_trim_microhomology(csv_tags, expected):
    result = trim_microhomology(csv_tags)
    assert result == expected, f"Expected {expected}, but got {result}"
