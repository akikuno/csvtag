from __future__ import annotations

import re


def split(csv_tag: str) -> list[str]:
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
    csv_tag = csv_tag.replace("cs:Z:", "")

    pattern = r"(\=[ACGTN]+|:[0-9]+|\*[ACGTN][ACGTN]|\+[ACGTN]+|\-[ACGTN]+|\~[ACGTN]{2}[0-9]+[ACGTN]{2}|[acgtn]+)"
    csv_split = [csv for csv in re.split(pattern, csv_tag) if csv]

    return csv_split
