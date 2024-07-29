from __future__ import annotations

import re
from typing import Generator


def split_by_tag(csv_tag: str) -> Generator[str]:
    """Split a csv tag
    Args:
        csv_tag (str): a csv tag
    Return:
        list[str]: splits a csv tag by operators

    Example:
        >>> import csvtag
        >>> csv_tag = ":4*AG:3"
        >>> csvtag.split(csv_tag)
        [':4', '*AG', ':3']
    """
    pattern = (
        r"(\=[ACGTN]+|:[0-9]+|\*[ACGTN][ACGTN]|\+[ACGTN]+|\-[ACGTN]+|\~[ACGTN]{2}[0-9]+[ACGTN]{2}|[=*+-~][acgtn]+)"
    )

    return (csv for csv in re.split(pattern, csv_tag) if csv)


def split_by_inversion(csv_tag: str) -> Generator[str]:
    """Split a csv tag by inversion
    Args:
        csv_tag (str): a csv tag
    Return:
        list[str]: splits a csv tag by inversion

    Example:
        >>> import csvtag
        >>> csv_tag = "=AA=aa*ga=a=AA"
        >>> csvtag.split_by_inversion(csv_tag)
        ['=AA', '=aa', '*ga', '=a', '=AA']
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
