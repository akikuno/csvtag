from __future__ import annotations

from collections.abc import Iterator
from itertools import groupby

import cstag

from csvtag.combiner import combine_splitted_csv_tag as combine
from csvtag.splitter import split_by_nucleotide as split

###########################################################
# remove_microholomogy
###########################################################


def _check_microhomology(curr_sequence: str, next_sequence: str) -> int:
    len_microhomology = 0
    min_length = min(len(curr_sequence), len(next_sequence))

    for i in range(1, min_length):
        if curr_sequence[-i:] == next_sequence[:i]:
            len_microhomology = i
        else:
            break
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

        # Check the length of microhomology
        len_microhomology = _check_microhomology(cstag.to_sequence(curr_cstag), cstag.to_sequence(next_cstag))
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
