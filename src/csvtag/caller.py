from __future__ import annotations

from collections.abc import Iterator
from itertools import groupby
from pathlib import Path

import cstag

from csvtag.microhomology_handler import remove_microhomology
from csvtag.sam_handler import extract_alignment, is_forward_strand, read_sam, remove_overlapped_alignments

############################
# On-Going!!!
############################

def convert_to_csvtag(alignments: list[dict[str, str | int]]) -> str:
    idx_inversion = set()
    idx = 0
    while idx + 1 < len(alignments):
        curr_align = alignments[idx]
        next_align = alignments[idx + 1]

        curr_pos: int = curr_align["POS"]
        next_pos: int = next_align["POS"]
        
        if curr_pos == next_pos:
            idx += 1
            continue
        
        curr_flag: int = curr_align["FLAG"]
        next_flag: int = next_align["FLAG"]
        curr_cstag: str = curr_align["CSTAG"]
        next_cstag: str = next_align["CSTAG"]

        if is_forward_strand(curr_flag) != is_forward_strand(next_flag):
            if is_forward_strand(curr_flag):
                curr_align["CSTAG"] = curr_cstag.upper()
                next_align["CSTAG"] = cstag.revcomp(next_cstag).lower()
                idx_inversion.add(idx + 1)
            else:
                curr_align["CSTAG"] = cstag.revcomp(curr_cstag).lower()
                next_align["CSTAG"] = next_cstag.upper()
                idx_inversion.add(idx)
            next_align["FLAG"] -= 16
        else:
            if idx not in idx_inversion:
                curr_align["CSTAG"] = curr_cstag.upper()
            if idx + 1 not in idx_inversion:
                next_align["CSTAG"] = next_cstag.upper()

        idx += 1

    return "".join(align["CSTAG"] for align in alignments)


def call_csvtag(path_sam: str | Path) -> Iterator[dict[str, list[str]]]:
    alignments: Iterator[dict[str, str | int]] = extract_alignment(read_sam(path_sam))
    alignments = remove_overlapped_alignments(alignments)
    alignments = remove_microhomology(alignments)

    alignments = list(alignments)
    alignments.sort(key=lambda x: (x["QNAME"], x["POS"]))
    for qname, alignments_grouped in groupby(alignments, key=lambda x: x["QNAME"]):
        alignments_grouped = list(alignments_grouped)
        if len(alignments_grouped) <= 2:
            yield {qname: [a["CSTAG"].upper() for a in alignments_grouped]}
            continue
        csv_tag = convert_to_csvtag(alignments_grouped)
        yield {qname: [csv_tag]}
