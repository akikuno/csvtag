from __future__ import annotations

from csvtag.splitter import split_by_nucleotide as split
from csvtag.to_sequence import to_sequence


def _get_length_of_microhomology(curr_sequence: str, next_sequence: str) -> int:
    len_microhomology = 0
    min_length = min(len(curr_sequence), len(next_sequence))

    for i in range(1, min_length + 1):
        if curr_sequence[-i:] == next_sequence[:i]:  # and curr_qual[-i:] == next_qual[:i]:
            len_microhomology = i
    return len_microhomology


###########################################################
# Handle csv tag
###########################################################


def trim_microhomology(csv_tags: list[str]) -> list[str]:
    """Trim microhomology from csv tags

    csv_tags = ["=AAATTT", "=TTTCCC"]
    trim_microhomology(csv_tags)
    ["=AAATTT", "=CCC"]
    """
    from csvtag.combiner import combine_splitted_tags

    csv_tags_trimmed = [csv_tags[0]]
    visited = {0}
    for i in range(len(csv_tags) - 1):
        curr_csvtag, next_csvtag = csv_tags[i], csv_tags[i + 1]

        # Check the length of microhomology
        len_microhomology = _get_length_of_microhomology(
            to_sequence(curr_csvtag.upper()), to_sequence(next_csvtag.upper())
        )

        csvtag_trimmed = list(split(next_csvtag))[len_microhomology:]
        next_csvtag = combine_splitted_tags(csvtag_trimmed)

        if i not in visited:
            csv_tags_trimmed.append(curr_csvtag)

        csv_tags_trimmed.append(next_csvtag)

        visited.add(i)
        visited.add(i + 1)

    return csv_tags_trimmed
