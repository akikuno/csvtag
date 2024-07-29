from __future__ import annotations

from csvtag.revcomp import revcomp
from csvtag.spliter import split_by_inversion, split_by_tag


def to_sequence(csv_tag: str) -> str:
    """Reconstruct the reference subsequence in the alignment

    Args:
        csv_tag (str): csv tag in the **long** format

    Returns:
        str: The sequence string derived from the cs tag.

    Example:
        >>> import csvtag
        >>> csv_tag = "=AA=aa*ga=a=AA"
        >>> csvtag.to_sequence(csv_tag)
        'AAttttAA'
    """
    sequence = []

    csv_tag_inversion_corrected = []
    for csv in split_by_inversion(csv_tag):
        nucleotide = csv[-1]
        if nucleotide.islower():
            csv_tag_inversion_corrected.append(revcomp(csv))
        else:
            csv_tag_inversion_corrected.append(csv)
    csv_tag_inversion_corrected = "".join(csv_tag_inversion_corrected)

    for cs in split_by_tag(csv_tag_inversion_corrected):
        if cs.startswith("="):
            sequence.append(cs[1:])
        elif cs.startswith("+"):
            sequence.append(cs[1:])
        elif cs.startswith("*"):
            sequence.append(cs[-1])

    return "".join(sequence)
