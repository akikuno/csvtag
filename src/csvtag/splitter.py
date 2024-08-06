from __future__ import annotations

import re
from collections.abc import Iterator


def split_by_tag(csv_tag: str) -> Iterator[str]:
    """Split a csv tag
    Args:
        csv_tag (str): a csv tag
    Return:
        list[str]: splits a csv tag by operators

    Example:
        >>> import csvtag
        >>> csv_tag = ":4*AG:3"
        >>> csvtag.split_by_tag(csv_tag)
        [':4', '*AG', ':3']
    """
    pattern = (
        r"(\=[ACGTN]+|:[0-9]+|\*[ACGTN][ACGTN]|\+[ACGTN]+|\-[ACGTN]+|\~[ACGTN]{2}[0-9]+[ACGTN]{2}|[=*+-~][acgtn]+)"
    )

    return (csv for csv in re.split(pattern, csv_tag) if csv)


def split_by_inversion(csv_tag: str) -> Iterator[str]:
    """Split a csv tag by inversion
    Args:
        csv_tag (str): a csv tag
    Return:
        list[str]: splits a csv tag by inversion

    Example:
        >>> import csvtag
        >>> csv_tag = "=AA=aa*ga=a=AA"
        >>> csvtag.split_by_inversion(csv_tag)
        ["=AA", "=aa*ga=a", "=TT"]
    """
    csv_tag_inversion = []
    is_inversion = False

    for csv in split_by_tag(csv_tag):
        nucleotide = csv[-1]
        if nucleotide.islower():
            csv_tag_inversion.append(csv)
            if not is_inversion:
                is_inversion = True

        elif is_inversion:
            yield ("".join(csv_tag_inversion))
            yield csv
            csv_tag_inversion = []
            is_inversion = False

        else:
            yield csv

    if csv_tag_inversion:  # if the last element is an inversion
        yield "".join(csv_tag_inversion)


###########################################################
# split_by_nucleotide
###########################################################


def handle_insertion(csv_tag: str, csv_tag_next: str) -> list[str]:
    """Handles insertion operations.
    csv_tag = "+acgt"
    csv_tag_next = "=AAA"
    results = handle_insertion(csv_tag, csv_tag_next)
    expected = ["+a|+c|+g|+t|=A", "=A", "=A"]
    """
    insertion: str = "|".join(handle_match_deletion(csv_tag, "+"))

    operand_next = csv_tag_next[0]
    if operand_next in {"=", "-"}:
        csv_tag_split = iter(handle_match_deletion(csv_tag_next, operand_next))
        return [insertion + "|" + next(csv_tag_split), *list(csv_tag_split)]
    elif operand_next == "*":
        return [insertion + "|" + csv_tag_next]
    else:
        n_bases = handle_splice(csv_tag_next)
        return [insertion + "|" + "=N"] + n_bases[1:]


def handle_splice(csv_tag: str) -> list[str]:
    """Handles splice operations."""
    match = re.match(r"([A-Za-z][A-Za-z])([0-9]+)([A-Za-z][A-Za-z])", csv_tag.replace("~", ""))
    if match is None:
        return []
    _, splice, _ = match.groups()
    return ["=N"] * int(splice)


def handle_match_deletion(csv_tag: str, operand: str) -> list[str]:
    """Handles substitution or deletion operations."""
    return [operand + c for c in csv_tag.replace(operand, "")]


def split_by_nucleotide(csv_tag: str) -> Iterator[str]:
    """Generate CS SPLIT, a comma-separated nucleotide sequence

    Args:
        csvtag (str): a long format csvtag

    Returns:
        str: csv split

    Examples:
        >>> csv_tag = "=A+TTT=CC-AA=T*AG=TT"
        >>> list(split_by_nucleotide(csv_tag))
        ["=A", "+T|+T|+T|=C", "=C", "-A", "-A", "=T", "*AG", "=T", "=T"

        >>> csv_tag = "=A~AA5CC=A"
        >>> list(split_by_nucleotide(csv_tag))
        ["=A", "=A", "=A", "=N", "=N", "=N", "=N", "=N", "=C", "=C", "=A"]

    """
    csv_tag_splitted = split_by_tag(csv_tag)
    csv_tags = []
    for csv_tag in csv_tag_splitted:
        operand = csv_tag[0]
        if operand == "*":
            csv_tags.append(csv_tag)
        elif operand == "+":
            csv_tags += handle_insertion(csv_tag, next(csv_tag_splitted))
        elif operand == "~":
            csv_tags += handle_splice(csv_tag)
        else:
            csv_tags += handle_match_deletion(csv_tag, operand)
    yield from csv_tags

