from __future__ import annotations

import re
from collections.abc import Iterator
from itertools import groupby
from pathlib import Path

import cstag


def read_sam(path_of_sam: str | Path) -> Iterator[list[str]]:
    with open(path_of_sam) as f:
        for line in f:
            yield line.strip().split("\t")


###########################################################
# Format headers and alignments
###########################################################


def extract_sqheaders(sam: list[list[str]]) -> dict[str, int]:
    """Extract SN (Reference sequence name) and LN (Reference sequence length) from SQ header

    Args:
        sam (list[list[str]]): a list of lists of SAM format

    Returns:
        dict: a dictionary containing (multiple) SN and LN
    """
    sqheaders = [s for s in sam if "@SQ" in s]
    sn_ln_output = {}
    for sqheader in sqheaders:
        sn_ln = [sq for sq in sqheader if re.search(("SN:|LN:"), sq)]
        sn = sn_ln[0].replace("SN:", "")
        ln = sn_ln[1].replace("LN:", "")
        sn_ln_output |= {sn: int(ln)}
    return sn_ln_output


def extract_alignment(sam: list[list[str]]) -> Iterator[dict[str, str | int]]:
    """Extract mapped alignments from SAM

    Args:
        sam (list[list[str]]): a list of lists of SAM format including cs tag

    Returns:
        Iterator[dict[str, str | int]]: a dictionary containing QNAME, FLAG, RNAME, POS, CIGAR, SEQ, QUAL, CSTAG
    """
    for alignment in sam:
        if alignment[0].startswith("@") or alignment[2] == "*" or alignment[9] == "*":
            continue
        idx_cstag = next((i for i, a in enumerate(alignment) if a.startswith("cs:Z:")), None)
        yield dict(
            QNAME=alignment[0].replace(",", "_"),
            FLAG=int(alignment[1]),
            RNAME=alignment[2],
            POS=int(alignment[3]),
            MAPQ=int(alignment[4]),
            CIGAR=alignment[5],
            SEQ=alignment[9],
            QUAL=alignment[10],
            CSTAG=alignment[idx_cstag].replace("cs:Z", ""),
        )


###########################################################
# Remove Softclips from SEQ and QUAL
###########################################################


def _split_cigar(cigar: str) -> Iterator[str]:
    cigar_iter = iter(re.split(r"([MIDNSHPX=])", cigar))
    return (i + op for i, op in zip(cigar_iter, cigar_iter))


def _calculate_alignment_length(cigar: str) -> int:
    return sum(int(c[:-1]) for c in _split_cigar(cigar) if c[-1] in "MDN=X")


def _is_complete_overlapped(prev_cigar: str, prev_pos: int, curr_cigar: str, curr_pos: int) -> bool:
    """Detect the shorter reads that are completely included in the longer reads"""
    prev_start = prev_pos - 1
    prev_end = prev_start + _calculate_alignment_length(prev_cigar)
    curr_start = curr_pos - 1
    curr_end = curr_start + _calculate_alignment_length(curr_cigar)

    if prev_start <= curr_start and prev_end >= curr_end:
        return True

    return False


def _is_non_microhomologic_overlapped(
    prev_cstag: str, curr_cstag, prev_cigar: str, prev_pos: int, curr_cigar: str, curr_pos: int
) -> bool:
    prev_seq = cstag.to_sequence(prev_cstag)
    curr_seq = cstag.to_sequence(curr_cstag)

    prev_start = prev_pos - 1
    prev_end = prev_start + _calculate_alignment_length(prev_cigar)
    curr_start = curr_pos - 1
    curr_end = curr_start + _calculate_alignment_length(curr_cigar)

    overlap_length = min(prev_end, curr_end) - max(prev_start, curr_start)

    for prev_base, curr_base in zip(prev_seq[::-1][:overlap_length], curr_seq[:overlap_length]):
        if prev_base != curr_base:
            return True

    return False


def _is_overlapped(prev_cstag: str, curr_cstag, prev_cigar: str, prev_pos: int, curr_cigar: str, curr_pos: int) -> bool:
    return _is_complete_overlapped(prev_cigar, prev_pos, curr_cigar, curr_pos) or _is_non_microhomologic_overlapped(
        prev_cstag, curr_cstag, prev_cigar, prev_pos, curr_cigar, curr_pos
    )


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
    alignments.sort(key=lambda x: [x["QNAME"], x["POS"]])
    alignments_group = groupby(alignments, lambda x: x["QNAME"])

    for _, alignments in alignments_group:
        alignments = list(alignments)
        if len(alignments) == 1:
            yield from (alignment for alignment in alignments)
            continue

        for i, (previous_alignment, current_alignment) in enumerate(zip(alignments, alignments[1:])):
            prev_cstag = previous_alignment["CSTAG"]
            curr_cstag = current_alignment["CSTAG"]
            prev_cigar = previous_alignment["CIGAR"]
            curr_cigar = current_alignment["CIGAR"]
            prev_pos = previous_alignment["POS"]
            curr_pos = current_alignment["POS"]

            if _is_overlapped(prev_cstag, curr_cstag, prev_cigar, prev_pos, curr_cigar, curr_pos):
                if _calculate_alignment_length(prev_cigar) >= _calculate_alignment_length(curr_cigar):
                    alignments[i + 1] = previous_alignment
                else:
                    alignments[i] = current_alignment

        yield from (alignment for alignment in _remove_duplicates(alignments))
