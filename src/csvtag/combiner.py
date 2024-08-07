from __future__ import annotations

from collections.abc import Iterator

from csvtag.microhomology_handler import trim_microhomology
from csvtag.splitter import split_by_tag
from csvtag.to_sequence import to_sequence

###########################################################
# combine_splitted_csv_tag
###########################################################


def _split_csv_tag_by_insertion(csv_tags: Iterator[str]) -> Iterator[str]:
    for csv_tag in csv_tags:
        if csv_tag.startswith("+"):  # Insertion
            yield from csv_tag.split("|")
        else:
            yield csv_tag


def combine_splitted_tags(splitted_csv_tag: Iterator[str]) -> str:
    """Conbine splitted csv tag

    Args:
        splitted csv tag (str)

    Returns:
        str: csv tag

    Examples:
        >>> splitted_csv_tag = iter(["=A", "+T|+T|+T|=C", "=C", "-A", "-A", "=T", "*AG", "=T", "=T"])
        >>> combine_splitted_csv_tag(splitted_csv_tag)
        "=A+TTT=CC-AA=T*AG=TT"

        >>> splitted_csv_tag = iter(["=A", "=A", "=A", "=N", "=N", "=N", "=N", "=N", "=C", "=C", "=A"])
        >>> combine_splitted_csv_tag(csv_tag)
        "=AAANNNNNCCA"
    """
    splitted_csv_tag = _split_csv_tag_by_insertion(splitted_csv_tag)

    combined_csv_tags = []
    prev_csv_tag = next(splitted_csv_tag)
    prev_prefix = prev_csv_tag[0]
    combined_csv = [prev_csv_tag[1:]]

    for curr_csv_tag in splitted_csv_tag:
        curr_prefix, curr_csv = curr_csv_tag[0], curr_csv_tag[1:]
        if prev_prefix == curr_prefix:
            if curr_prefix == "*":
                combined_csv.append(curr_csv_tag)
            else:
                combined_csv.append(curr_csv)
        else:
            combined_csv_tags.append(prev_prefix + "".join(combined_csv))
            prev_prefix = curr_prefix
            combined_csv = [curr_csv]

    combined_csv_tags.append(prev_prefix + "".join(combined_csv))
    return "".join(combined_csv_tags)


###########################################################
# combine_neighboring_csv_tags
###########################################################


def _get_n_lengths(csv_tags: list[str], positions: list[int]) -> list[int]:
    n_lengths = []
    for curr_tag, _, curr_pos, next_pos in zip(csv_tags, csv_tags[1:], positions, positions[1:]):
        n_length = next_pos - curr_pos - len(to_sequence(curr_tag))
        n_lengths.append(n_length)
    n_lengths.append(-1)
    return n_lengths


def _group_tags_by_distance(csv_tags: list[str], n_lengths: list[int], distance: int) -> list[list[str]]:
    groups = []
    current_group = [csv_tags[0]]

    for i in range(len(n_lengths) - 1):
        if n_lengths[i] <= distance:
            current_group.append(csv_tags[i + 1])
        else:
            groups.append(current_group)
            current_group = [csv_tags[i + 1]]

    groups.append(current_group)
    return groups


def _padding_n(csv_tag: str, n_length: int, side: str = "left") -> str:
    csv_tag_split = list(split_by_tag(csv_tag))
    is_inversion = csv_tag_split[0].islower()
    n_character = "n" if is_inversion else "N"

    if side == "left":
        if csv_tag_split[0].startswith("="):
            return "=" + (n_character * n_length) + csv_tag.lstrip("=")
        else:
            return "=" + (n_character * n_length) + csv_tag
    else:
        if csv_tag_split[-1].startswith("="):
            return csv_tag + (n_character * n_length)
        else:
            return csv_tag + "=" + (n_character * n_length)


def _combine_group(csv_tags: list[str], n_lengths: list[int], distance: int) -> list[str]:
    idx = 0
    tags_n_appended = []
    for group_csvtags in _group_tags_by_distance(csv_tags, n_lengths, distance):
        trimmed_csvtags = trim_microhomology(group_csvtags)
        n_length = n_lengths[idx : idx + len(trimmed_csvtags) - 1] + [-1]
        n_appended = []
        for tag, n in zip(trimmed_csvtags, n_length):
            if n != -1:
                n_appended.append(_padding_n(tag, n, side="right"))
            else:
                n_appended.append(tag)
        tags_n_appended.append(("".join(n_appended)))

        idx += len(trimmed_csvtags)

    return tags_n_appended


def combine_neighboring_csv_tags(csv_tags: list[str], positions: list[int], distance: int = 50) -> list[str]:
    """
    Examples:
        >>> csv_tags = ["=AA", "=tt", "=CC"]
        >>> positions = [1, 5, 9]
        >>> combine_neighboring_csv_tags(csv_tags, positions, 50)
        ['=AANN=tt=NNCC']

        >>> csv_tags = ["=AA", "=tt", "=CC"]
        >>> positions = [1, 5, 100]
        >>> combine_neighboring_csv_tags(csv_tags, positions, 50)
        ['=AANN=tt', '=CC']
    """
    if not csv_tags or not positions or len(csv_tags) != len(positions):
        raise ValueError("csv_tags and positions must be non-empty and of the same length.")

    n_lengths = _get_n_lengths(csv_tags, positions)

    return _combine_group(csv_tags, n_lengths, distance)
