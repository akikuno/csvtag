from __future__ import annotations

from collections.abc import Iterator
from itertools import groupby
from pathlib import Path

import cstag

from csvtag.microhomology_handler import remove_microhomology
from csvtag.sam_handler import (
    calculate_alignment_length,
    extract_alignment,
    is_forward_strand,
    read_sam,
    remove_overlapped_alignments,
)


def _is_second_strand_different(first_flag: int, second_flag: int, third_flag: int) -> bool:
    if is_forward_strand(first_flag) == is_forward_strand(third_flag) and is_forward_strand(
        first_flag
    ) != is_forward_strand(second_flag):
        return True
    else:
        return False


def _is_within_bases(first_end: int, second_pos: int, second_end: int, third_pos: int, base_num: int = 100) -> bool:
    if second_pos - first_end <= base_num and third_pos - second_end <= base_num:
        return True
    else:
        return False


def _padding_n(cs_tag: str, length: int, side: str = "left") -> str:
    if side == "left":
        if cstag.split(cs_tag)[0].startswith("="):
            return "=" + ("N" * length) + cs_tag.lstrip("=")
        else:
            return "=" + ("N" * length) + cs_tag
    else:
        if cstag.split(cs_tag)[-1].startswith("="):
            return cs_tag + ("N" * length)
        else:
            return cs_tag + "=" + ("N" * length)


def _unique_dicts(list_of_dicts: list[dict[str, str]]) -> list[dict[str, str]]:
    unique_items = {tuple(d.items()) for d in list_of_dicts}
    unique_dicts = [dict(t) for t in unique_items]
    return unique_dicts


def convert_to_csvtag(alignments: list[dict[str, str | int]]) -> list[dict[str, str | int]]:
    idx = 0
    while idx + 2 < len(alignments):
        first_align = alignments[idx]
        second_align = alignments[idx + 1]
        third_align = alignments[idx + 2]

        # (1) When only the Flag of the second among first, second, and third is different
        # (2) When the distance between the End of first and the Start of second, and the distance between the End of second and the Start of third are both within certain bases
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

        if is_second_strand_different and is_within_bases:
            first_cstag: str = first_align["CSTAG"]
            second_cstag: str = second_align["CSTAG"]
            third_cstag: str = third_align["CSTAG"]

            if second_pos - first_end > 0:
                first_cstag = _padding_n(first_cstag, second_pos - first_end, side="right")
            if third_pos - second_end > 0:
                third_cstag = _padding_n(third_cstag, third_pos - second_end, side="left")

            csv_tag = first_cstag.upper() + second_cstag.lower() + third_cstag.upper()

            alignments[idx]["POS"] = first_pos
            alignments[idx + 1]["POS"] = first_pos
            alignments[idx + 2]["POS"] = first_pos

            alignments[idx]["CSTAG"] = csv_tag
            alignments[idx + 1]["CSTAG"] = csv_tag
            alignments[idx + 2]["CSTAG"] = csv_tag

        idx += 1

    return alignments


def _revcomp_cstag_of_reverse_strand(alignments: list[dict[str, str]]) -> list[dict[str, str]]:
    for alignment in alignments:
        if not is_forward_strand(alignment["FLAG"]):
            alignment["CSTAG"] = cstag.revcomp(alignment["CSTAG"])
    return alignments


def _upper_cstag(alignments: list[dict[str, str]]) -> list[dict[str, str]]:
    for alignment in alignments:
        alignment["CSTAG"] = alignment["CSTAG"].upper()
    return alignments


def call_csvtag(path_sam: str | Path) -> Iterator[dict[str, str | int]]:
    """
    output: [{"QNAME": "read1", "CSVTAG": "=AAAAA"], {"QNAME": "read1", "CSVTAG": "=TTTTT"}, ...]}
    """
    alignments: Iterator[dict[str, str | int]] = extract_alignment(read_sam(path_sam))
    alignments = remove_overlapped_alignments(alignments)
    alignments = remove_microhomology(alignments)

    alignments = list(alignments)
    alignments.sort(key=lambda x: (x["QNAME"], x["RNAME"], x["POS"]))
    for (qname, rname), alignments_grouped in groupby(alignments, key=lambda x: [x["QNAME"], x["RNAME"]]):
        alignments_grouped = list(alignments_grouped)

        # Convert all cs tags to the plus strand
        alignments_grouped = _revcomp_cstag_of_reverse_strand(alignments_grouped)

        # Convert all CSV tags to uppercase (Note: they will no longer be standard cs tags)
        alignments_grouped = _upper_cstag(alignments_grouped)

        if len(alignments_grouped) <= 2:
            yield from [
                {"QNAME": qname, "RNAME": rname, "POS": a["POS"], "CSVTAG": a["CSTAG"]} for a in alignments_grouped
            ]
            continue

        alignments_grouped = convert_to_csvtag(alignments_grouped)

        yield from _unique_dicts(
            [{"QNAME": qname, "RNAME": rname, "POS": a["POS"], "CSVTAG": a["CSTAG"]} for a in alignments_grouped]
        )
