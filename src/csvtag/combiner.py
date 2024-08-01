from __future__ import annotations

from collections.abc import Iterator


def _split_csv_tag_by_insertion(csv_tag: Iterator[str]) -> Iterator[str]:
    for csv in csv_tag:
        if csv.startswith("+"):  # Insertion
            yield from csv.split("|")
        else:
            yield csv


def combine_splitted_csv_tag(splitted_csv_tag: Iterator[str]) -> str:
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



