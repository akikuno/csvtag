from __future__ import annotations

from dataclasses import dataclass

import pytest

from csvtag.overlap_handler import (
    _is_complete_overlapped,
    # _is_non_microhomologic_overlapped,
    remove_overlapped_alignments,
)

#########################################
# _is_overlapped
#########################################


@dataclass
class OverlappedAlignment:
    curr_cigar: str
    next_cigar: str
    curr_cstag: str = ""
    next_cstag: str = ""
    curr_pos: int = 0
    next_pos: int = 0


@pytest.mark.parametrize(
    "alignments_overlapped, expected",
    [
        (OverlappedAlignment("10M", "5M", curr_pos=1, next_pos=1), True),
        (OverlappedAlignment("10M", "5M", curr_pos=1, next_pos=6), True),
        (OverlappedAlignment("10M", "5M", curr_pos=1, next_pos=7), False),
        (OverlappedAlignment("5M", "10M", curr_pos=1, next_pos=1), False),
        (OverlappedAlignment("10M", "10M", curr_pos=1, next_pos=1), True),
        (OverlappedAlignment("10M5I10M", "5M", curr_pos=1, next_pos=11), True),
        (OverlappedAlignment("10M5I10M", "15M", curr_pos=1, next_pos=11), False),
    ],
)
def test_is_complete_overlapped(alignments_overlapped, expected):
    result = _is_complete_overlapped(alignments_overlapped)
    assert result == expected, f"Expected {expected}, but got {result}"


# @pytest.mark.parametrize(
#     "alignments_overlapped, expected",
#     [
#         pytest.param(
#             OverlappedAlignment(
#                 curr_cstag="=GACT", next_cstag="=ACT", curr_cigar="4M", next_cigar="3M", curr_pos=1, next_pos=1
#             ),
#             True,
#             id="test_case_1",
#         ),
#         pytest.param(
#             OverlappedAlignment(
#                 curr_cstag="=ACTG", next_cstag="=ACT", curr_cigar="4M", next_cigar="3M", curr_pos=1, next_pos=6
#             ),
#             False,
#             id="test_case_2",
#         ),
#         pytest.param(
#             OverlappedAlignment(
#                 curr_cstag="=ACTG", next_cstag="=ACT", curr_cigar="4M", next_cigar="3M", curr_pos=1, next_pos=7
#             ),
#             False,
#             id="test_case_3",
#         ),
#         pytest.param(
#             OverlappedAlignment(
#                 curr_cstag="=ACTG", next_cstag="=ACTG", curr_cigar="4M", next_cigar="4M", curr_pos=1, next_pos=1
#             ),
#             True,
#             id="test_case_5",
#         ),
#         pytest.param(
#             OverlappedAlignment(
#                 curr_cstag="=ACTGACTGACTG",
#                 next_cstag="=ACTGACTGACTG",
#                 curr_cigar="12M",
#                 next_cigar="12M",
#                 curr_pos=1,
#                 next_pos=11,
#             ),
#             False,
#             id="test_case_6",
#         ),
#     ],
# )
# def test_is_non_microhomologic_overlapped(alignments_overlapped, expected):
#     result = _is_non_microhomologic_overlapped(alignments_overlapped)
#     assert result == expected, f"Expected {expected}, but got {result}"


@pytest.mark.parametrize(
    "alignments, expected",
    [
        pytest.param(
            [
                {"RNAME": "chr1", "QNAME": "read1", "POS": 1, "CIGAR": "10M", "CSTAG": "=ACTGACTGAC"},
                {"RNAME": "chr1", "QNAME": "read2", "POS": 1, "CIGAR": "10M", "CSTAG": "=TGCATGCATG"},
            ],
            [
                {"RNAME": "chr1", "QNAME": "read1", "POS": 1, "CIGAR": "10M", "CSTAG": "=ACTGACTGAC"},
                {"RNAME": "chr1", "QNAME": "read2", "POS": 1, "CIGAR": "10M", "CSTAG": "=TGCATGCATG"},
            ],
            id="no_duplicates",
        ),
        pytest.param(
            [
                {"RNAME": "chr1", "QNAME": "read1", "POS": 1, "CIGAR": "10M", "CSTAG": "=ACTGACTGAC"},
                {"RNAME": "chr1", "QNAME": "read1", "POS": 1, "CIGAR": "5M", "CSTAG": "=ACTGA"},
            ],
            [{"RNAME": "chr1", "QNAME": "read1", "POS": 1, "CIGAR": "10M", "CSTAG": "=ACTGACTGAC"}],
            id="completely_overlapped_reads",
        ),
        pytest.param(
            [
                {"RNAME": "chr1", "QNAME": "read1", "POS": 1, "CIGAR": "5M", "CSTAG": "=ACTGA"},
                {"RNAME": "chr1", "QNAME": "read1", "POS": 3, "CIGAR": "5M", "CSTAG": "=GACCT"},
            ],
            [
                {"RNAME": "chr1", "QNAME": "read1", "POS": 1, "CIGAR": "5M", "CSTAG": "=ACTGA"},
                {"RNAME": "chr1", "QNAME": "read1", "POS": 3, "CIGAR": "5M", "CSTAG": "=GACCT"},
            ],
            id="microhomology_reads",
        ),
    ],
)
def test_remove_overlapped_alignments(alignments, expected):
    result = list(remove_overlapped_alignments(iter(alignments)))
    result = sorted(result, key=lambda x: [x["QNAME"], x["RNAME"], x["POS"]])
    assert result == expected, f"Expected {expected}, but got {result}"
