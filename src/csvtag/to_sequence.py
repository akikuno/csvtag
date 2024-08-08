from __future__ import annotations

from csvtag.splitter import split_by_tag


def to_sequence(csv_tag: str) -> str:
    """Reconstruct the **query** subsequence in the alignment

    Args:
        csv_tag (str): csv tag

    Returns:
        str: The sequence string derived from the csv tag.

    Example:
        >>> import csvtag
        >>> csv_tag = "=AA=aa*ga=A-CC+G|G|G|=AAA"
        >>> csvtag.to_sequence(csv_tag)
        'AAaaaAGGGAAAA'
    """
    sequence = []

    for tag in split_by_tag(csv_tag):
        if tag.startswith("="):
            sequence.append(tag[1:])
        elif tag.startswith("+"):
            sequence.append(tag[1:])
        elif tag.startswith("*"):
            sequence.append(tag[-1])

    return "".join(sequence)
