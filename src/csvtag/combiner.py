from __future__ import annotations

from collections.abc import Iterator

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


def _padding_n(csv_tag: str, n_length: int, side: str = "left") -> str:
    csv_tag_split = list(split_by_tag(csv_tag))
    if side == "left":
        if csv_tag_split[0].startswith("="):
            return "=" + ("N" * n_length) + csv_tag.lstrip("=")
        else:
            return "=" + ("N" * n_length) + csv_tag
    else:
        if csv_tag_split[-1].startswith("="):
            return csv_tag + ("N" * n_length)
        else:
            return csv_tag + "=" + ("N" * n_length)


def combine_neighboring_csv_tags(csv_tags: list[str], positions: list[int], distance: int = 50) -> Iterator[str]:
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

    # visited = set()
    i = 0
    combined_csv_tags = []
    while i + 1 < len(csv_tags):
        curr_csv_tag = csv_tags[i]
        next_csv_tag = csv_tags[i + 1]
        curr_pos = positions[i]
        next_pos = positions[i + 1]

        n_length = next_pos - curr_pos - len(to_sequence(curr_csv_tag))

        if 0 < n_length < distance:
            combined_csv_tags.append(_padding_n(curr_csv_tag, n_length, side="right"))
            # visited.add(i)
            # visited.add(i + 1)
        else:
            if combined_csv_tags:
                yield "".join(combined_csv_tags)
                combined_csv_tags = []
            # if i not in visited:
            # yield curr_csv_tag
        i += 1

    # if i + 1 not in visited:
    if combined_csv_tags:
        yield "".join(combined_csv_tags)
    for next_csv_tag in csv_tags[i:]:
        yield next_csv_tag
