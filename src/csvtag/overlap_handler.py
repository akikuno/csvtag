from collections.abc import Iterator
from dataclasses import dataclass
from itertools import groupby

from csvtag.sam_handler import calculate_alignment_length

###########################################################
# Remove Overlapped alignments
###########################################################


@dataclass
class OverlappedAlignment:
    curr_cigar: str
    next_cigar: str
    curr_cstag: str
    next_cstag: str
    curr_pos: int
    next_pos: int


def _is_complete_overlapped(alignments_overlapped: OverlappedAlignment) -> bool:
    """Detect the shorter reads that are completely included in the longer reads"""
    curr_pos = alignments_overlapped.curr_pos
    next_pos = alignments_overlapped.next_pos
    curr_cigar = alignments_overlapped.curr_cigar
    next_cigar = alignments_overlapped.next_cigar

    curr_start = curr_pos - 1
    curr_end = curr_start + calculate_alignment_length(curr_cigar)
    next_start = next_pos - 1
    next_end = next_start + calculate_alignment_length(next_cigar)

    if curr_start <= next_start and curr_end >= next_end:
        return True

    return False


# def _is_non_microhomologic_overlapped(alignments_overlapped: OverlappedAlignment) -> bool:
#     curr_pos = alignments_overlapped.curr_pos
#     next_pos = alignments_overlapped.next_pos
#     curr_cigar = alignments_overlapped.curr_cigar
#     next_cigar = alignments_overlapped.next_cigar
#     curr_cstag = alignments_overlapped.curr_cstag
#     next_cstag = alignments_overlapped.next_cstag

#     curr_seq = cstag.to_sequence(curr_cstag)
#     next_seq = cstag.to_sequence(next_cstag)

#     curr_start = curr_pos - 1
#     curr_end = curr_start + calculate_alignment_length(curr_cigar)
#     next_start = next_pos - 1
#     next_end = next_start + calculate_alignment_length(next_cigar)

#     overlap_length = min(curr_end, next_end) - max(curr_start, next_start)

#     if overlap_length <= 0:
#         return False

#     if curr_seq[-overlap_length:] == next_seq[:overlap_length]:
#         return False
#     else:
#         True


def _is_overlapped(alignments_overlapped: OverlappedAlignment) -> bool:
    return _is_complete_overlapped(alignments_overlapped)  # or _is_non_microhomologic_overlapped(alignments_overlapped)


def _remove_duplicates(list_of_dicts: list[dict[str, str | int]]) -> list[dict[str, str | int]]:
    """
    Remove duplicate dictionaries from a list of dictionaries.

    Args:
        list_of_dicts (list[dict[str, Any]]): A list of dictionaries.

    Returns:
        list[dict[str, Any]]: A list of dictionaries with duplicates removed.
    """
    unique_dicts = list({tuple(sorted(d.items())) for d in list_of_dicts})
    return [dict(t) for t in unique_dicts]


def remove_overlapped_alignments(alignments: Iterator[dict[str, str | int]]) -> Iterator[dict[str, str | int]]:
    """Remove non-microhomologic overlapped reads within the same QNAME.
    The overlapped sequences can be (1) realignments by microhomology or (2) resequence by sequencing error.
    The 'realignments' is not sequencing errors, and it preserves the same sequence.
    In contrast, the 'resequence' is a sequencing error with the following characteristics:
    (1) The shorter reads that are completely included in the longer reads
    (2) Overlapped but not the same DNA sequence
    The resequenced fragments will be discarded and the longest alignment will be retain.
    Example reads are in `tests/data/overlap/real_overlap.sam` and `tests/data/overlap/real_overlap2.sam`

    Args:
        alignments (list[dict[str, str | int]]): disctionalized alignments

    Returns:
        Iterator[dict[str, str | int]]: disctionalized alignments without overlaped reads
    """
    alignments = list(alignments)

    alignments.sort(key=lambda x: [x["QNAME"], x["RNAME"], x["POS"]])
    for _, alignments in groupby(alignments, lambda x: [x["QNAME"], x["RNAME"]]):
        alignments = list(alignments)
        if len(alignments) == 1:
            yield from (alignment for alignment in alignments)
            continue

        for i, (current_alignment, next_alignment) in enumerate(zip(alignments, alignments[1:])):
            alignments_overlapped = OverlappedAlignment(
                curr_cigar=current_alignment["CIGAR"],
                next_cigar=next_alignment["CIGAR"],
                curr_cstag=current_alignment["CSTAG"],
                next_cstag=next_alignment["CSTAG"],
                curr_pos=current_alignment["POS"],
                next_pos=next_alignment["POS"],
            )

            if _is_overlapped(alignments_overlapped):
                curr_length = calculate_alignment_length(alignments_overlapped.curr_cigar)
                next_length = calculate_alignment_length(alignments_overlapped.next_cigar)
                if curr_length >= next_length:
                    alignments[i + 1] = current_alignment
                else:
                    alignments[i] = next_alignment

        yield from (alignment for alignment in _remove_duplicates(alignments))
