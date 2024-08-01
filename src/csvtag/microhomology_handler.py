from __future__ import annotations

from collections.abc import Iterator
from itertools import groupby

import cstag

from csvtag.combiner import combine_splitted_csv_tag as combine
from csvtag.sam_handler import split_cigar
from csvtag.splitter import split_by_nucleotide as split

###########################################################
# remove_microholomogy
###########################################################


def _trim_softclip(qual: str, cigar: str) -> str:
    cigar_split = list(split_cigar(cigar))
    if cigar_split[0].endswith("S"):
        softclip_length = int(cigar_split[0][:-1])
        qual = qual[softclip_length:]
    if cigar_split[-1].endswith("S"):
        softclip_length = int(cigar_split[-1][:-1])
        qual = qual[:-softclip_length]
    return qual


def _get_length_of_microhomology(curr_sequence: str, next_sequence: str, curr_qual: str, next_qual: str) -> int:
    len_microhomology = 0
    min_length = min(len(curr_sequence), len(next_sequence))

    for i in range(1, min_length + 1):
        if curr_sequence[-i:] == next_sequence[:i] and curr_qual[-i:] == next_qual[:i]:
            len_microhomology = i
    return len_microhomology


def _trim_microhomology(alignments: list[dict[str, str | int]]) -> list[dict[str, str | int]]:
    """
    alignments = [{"QNAME": read1, "CSTAG": "=AAATTT"}, {"QNAME": read1, "CSTAG": "=TTTCCC"}]
    _trim_microhomology(alignments)
    [{"CSTAG": "=AAATTT"}, {"CSTAG": "=CCC"}]
    """
    idx = 0
    while idx < len(alignments) - 1:
        curr_align, next_align = alignments[idx], alignments[idx + 1]
        curr_cstag: str = curr_align["CSTAG"]
        next_cstag: str = next_align["CSTAG"]

        curr_cigar: str = curr_align["CIGAR"]
        next_cigar: str = next_align["CIGAR"]
        curr_qual: str = curr_align["QUAL"]
        next_qual: str = next_align["QUAL"]

        curr_qual = _trim_softclip(curr_qual, curr_cigar)
        next_qual = _trim_softclip(next_qual, next_cigar)

        # Check the length of microhomology
        len_microhomology = _get_length_of_microhomology(
            cstag.to_sequence(curr_cstag), cstag.to_sequence(next_cstag), curr_qual, next_qual
        )
        # Continue if no microhomology exists
        if len_microhomology == 0:
            idx += 1
            continue
        # Use less mutated alignment
        curr_microhomology = list(split(curr_cstag))[-len_microhomology:]
        next_microhomology = list(split(next_cstag))[:len_microhomology]
        matches_in_curr = sum(True if cs.startswith("=") and cs != "=N" else 0 for cs in curr_microhomology)
        matches_in_next = sum(True if cs.startswith("=") and cs != "=N" else 0 for cs in next_microhomology)
        # Remove frequently mutated microhomology region from the alignment's CSTAG
        if matches_in_curr >= matches_in_next:
            cstag_trimmed = list(split(next_cstag))[len_microhomology:]
            alignments[idx + 1]["CSTAG"] = combine(cstag_trimmed)
        else:
            cstag_trimmed = list(split(curr_cstag))[:-len_microhomology]
            alignments[idx]["CSTAG"] = combine(cstag_trimmed)
        idx += 1

    return alignments


def remove_microhomology(alignments: Iterator[dict[str, str | int]]) -> Iterator[dict[str, str | int]]:
    alignments = list(alignments)
    alignments.sort(key=lambda x: (x["QNAME"], x["POS"]))

    for _, alignments_grouped in groupby(alignments, key=lambda x: x["QNAME"]):
        alignments_grouped = list(alignments_grouped)
        if len(alignments_grouped) == 1:
            yield from alignments_grouped
            continue
        yield from _trim_microhomology(alignments_grouped)
