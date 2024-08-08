from __future__ import annotations

from collections.abc import Iterator
from itertools import groupby
from pathlib import Path

import cstag

from csvtag.overlap_remover import remove_overlapped_alignments
from csvtag.sam_handler import (
    calculate_alignment_length,
    extract_alignment,
    is_forward_strand,
    read_sam,
)


def _is_second_strand_different(first_flag: int, second_flag: int, third_flag: int) -> bool:
    if is_forward_strand(first_flag) == is_forward_strand(third_flag) and is_forward_strand(
        first_flag
    ) != is_forward_strand(second_flag):
        return True
    else:
        return False


def _is_within_bases(
    first_end: int,
    second_pos: int,
    second_end: int,
    third_pos: int,
    base_num: int = 100,
) -> bool:
    if second_pos - first_end <= base_num and third_pos - second_end <= base_num:
        return True
    else:
        return False


def convert_to_csvtag(
    alignments: list[dict[str, str | int]],
) -> Iterator[dict[str, str | int]]:
    idx = 0
    visited = set()
    while idx + 2 < len(alignments):
        first_align = alignments[idx]
        second_align = alignments[idx + 1]
        third_align = alignments[idx + 2]

        # (1) When only the Flag of the second among first, second, and third is different
        # (2) When the distance between the End of first and the Start of second,
        #     and the distance between the End of second and the Start of third are both within certain bases
        # If both (1) and (2) are true, it is considered an inversion. Insert =N for regions that are apart.

        first_flag: int = first_align["FLAG"]
        second_flag: int = second_align["FLAG"]
        third_flag: int = third_align["FLAG"]

        first_pos: int = first_align["POS"]
        second_pos: int = second_align["POS"]
        third_pos: int = third_align["POS"]

        first_cigar: str = first_align["CIGAR"]
        second_cigar: str = second_align["CIGAR"]

        first_end: int = first_pos + calculate_alignment_length(first_cigar)
        second_end: int = second_pos + calculate_alignment_length(second_cigar)

        is_second_strand_different = _is_second_strand_different(first_flag, second_flag, third_flag)

        is_within_bases = _is_within_bases(first_end, second_pos, second_end, third_pos, base_num=50)

        first_cstag: str = first_align["CSTAG"]
        second_cstag: str = second_align["CSTAG"]

        if is_second_strand_different and is_within_bases:
            yield {
                "QNAME": first_align["QNAME"],
                "RNAME": first_align["RNAME"],
                "POS": first_pos,
                "CSVTAG": first_cstag,
            }
            yield {
                "QNAME": second_align["QNAME"],
                "RNAME": second_align["RNAME"],
                "POS": second_pos,
                "CSVTAG": second_cstag.lower(),
            }

            visited.add(idx)
            visited.add(idx + 1)
            idx += 2

        else:
            if idx not in visited:
                yield {
                    "QNAME": first_align["QNAME"],
                    "RNAME": first_align["RNAME"],
                    "POS": first_pos,
                    "CSVTAG": first_cstag,
                }

            visited.add(idx)
            idx += 1

    for i in range(idx, len(alignments)):
        alignment = alignments[i]
        yield {
            "QNAME": alignment["QNAME"],
            "RNAME": alignment["RNAME"],
            "POS": alignment["POS"],
            "CSVTAG": alignment["CSTAG"],
        }


def _revcomp_cstag_of_reverse_strand(
    alignments: list[dict[str, str | int]],
) -> list[dict[str, str]]:
    for alignment in alignments:
        if not is_forward_strand(alignment["FLAG"]):
            alignment["CSTAG"] = cstag.revcomp(alignment["CSTAG"])
    return alignments


def _upper_cstag(alignments: list[dict[str, str | int]]) -> list[dict[str, str | int]]:
    for alignment in alignments:
        alignment["CSTAG"] = alignment["CSTAG"].upper()
    return alignments


def call_csvtag(path_sam: str | Path) -> Iterator[dict[str, str | int]]:
    """
    Process SAM file and yield alignment information with CSV tags.

    This function reads a SAM file, extracts alignments, removes overlapped alignments,
    sorts them, and processes the alignment tags. The result is an iterator of dictionaries,
    each containing the query name, reference name, position, and processed CSV tag.

    Args:
        path_sam (str | Path): The path to the SAM file to be processed.

    Yields:
        Iterator[dict[str, str | int]]: An iterator of dictionaries with the following keys:
            - "QNAME" (str): Query name (read name).
            - "RNAME" (str): Reference sequence name.
            - "POS" (int): 1-based leftmost mapping position.
            - "CSVTAG" (str): Processed CSV tag in uppercase.

    Example:
        >>> for alignment in call_csvtag("example.sam"):
        ...     print(alignment)
        {"QNAME": "read1", "RNAME": "chr1", "POS": 100, "CSVTAG": "=AAAAA"}
        {"QNAME": "read1", "RNAME": "chr1", "POS": 150, "CSVTAG": "=TTTTT"}
        ...
    """
    alignments: Iterator[dict[str, str | int]] = extract_alignment(read_sam(path_sam))
    alignments = remove_overlapped_alignments(alignments)

    alignments = list(alignments)
    alignments.sort(key=lambda x: (x["QNAME"], x["RNAME"], x["POS"]))
    for (qname, rname), alignments_grouped in groupby(alignments, key=lambda x: [x["QNAME"], x["RNAME"]]):
        alignments_grouped = list(alignments_grouped)

        # Convert all cs tags to the plus strand
        alignments_grouped = _revcomp_cstag_of_reverse_strand(alignments_grouped)

        # Convert all CSV tags to uppercase (Note: they will no longer be standard cs tags)
        alignments_grouped = _upper_cstag(alignments_grouped)

        if len(alignments_grouped) <= 2:
            for alignment in alignments_grouped:
                yield {
                    "QNAME": qname,
                    "RNAME": rname,
                    "POS": alignment["POS"],
                    "CSVTAG": alignment["CSTAG"],
                }
            continue

        yield from convert_to_csvtag(alignments_grouped)
