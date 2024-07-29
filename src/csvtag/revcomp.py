from __future__ import annotations

import re

map_revcomp = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "n": "n",
}


def _extract_numbers(strings: str) -> list[str]:
    # Using regular expression to find all numbers in the strings: st
    return re.findall(r"\d+", strings)


def revcomp(csv_tag: str) -> str:
    """Converts a csv tag into its reverse complement.
    Args:
        csv_tag (str): a csv tag
        prefix (bool, optional): Whether to add the prefix 'cs:Z:' to the csv tag. Defaults to False

    Return:
        str: reverse complement of a csv tag

    Example:
        >>> import csvtag
        >>> csv_tag = "=AA=aa*ga=a=AA"
        >>> csvtag.revcomp(csv_tag)
        '=TT=t*ct=tt=TT'
    """
    pattern = (
        r"(\=[ACGTN]+|:[0-9]+|\*[ACGTN][ACGTN]|\+[ACGTN]+|\-[ACGTN]+|\~[ACGTN]{2}[0-9]+[ACGTN]{2}|[=*+-~][acgtn]+)"
    )

    csv_tag_revcomp = []
    for csv in re.split(pattern, csv_tag)[::-1]:
        if csv == "":
            continue
        operand = csv[0]
        if operand == ":":
            csv_tag_revcomp.append(csv)
        elif operand == "*":
            csv_tag_revcomp.append(f"*{map_revcomp[csv[1]]}{map_revcomp[csv[2]]}")
        elif operand == "~":
            numbers = _extract_numbers(csv)
            csv_tag_revcomp.append(
                f"~{map_revcomp[csv[-1]]}{map_revcomp[csv[-2]]}{numbers[0]}{map_revcomp[csv[2]]}{map_revcomp[csv[1]]}"
            )
        else:
            cs_revcomp = "".join([map_revcomp[c] for c in csv[1:]])[::-1]
            csv_tag_revcomp.append(f"{operand}{cs_revcomp}")

    return "".join(csv_tag_revcomp)
