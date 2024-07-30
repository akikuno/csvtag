from __future__ import annotations

import re
from collections.abc import Iterator

###########################################################
# Format headers and alignments
###########################################################


def extract_sqheaders(sam: list[list[str]]) -> dict[str, int]:
    """Extract SN (Reference sequence name) and LN (Reference sequence length) from SQ header

    Args:
        sam (list[list[str]]): a list of lists of SAM format

    Returns:
        dict: a dictionary containing (multiple) SN and LN
    """
    sqheaders = [s for s in sam if "@SQ" in s]
    sn_ln_output = {}
    for sqheader in sqheaders:
        sn_ln = [sq for sq in sqheader if re.search(("SN:|LN:"), sq)]
        sn = sn_ln[0].replace("SN:", "")
        ln = sn_ln[1].replace("LN:", "")
        sn_ln_output |= {sn: int(ln)}
    return sn_ln_output


def extract_alignment(sam: list[list[str]]) -> Iterator[dict[str, str | int]]:
    """Extract mapped alignments from SAM

    Args:
        sam (list[list[str]]): a list of lists of SAM format including cs tag

    Returns:
        Iterator[dict[str, str | int]]: a dictionary containing QNAME, FLAG, RNAME, POS, CIGAR, SEQ, QUAL, CSTAG
    """
    for alignment in sam:
        if alignment[0].startswith("@") or alignment[2] == "*" or alignment[9] == "*":
            continue
        idx_cstag = next((i for i, a in enumerate(alignment) if a.startswith("cs:Z:")), None)
        yield dict(
            QNAME=alignment[0].replace(",", "_"),
            FLAG=int(alignment[1]),
            RNAME=alignment[2],
            POS=int(alignment[3]),
            MAPQ=int(alignment[4]),
            CIGAR=alignment[5],
            SEQ=alignment[9],
            QUAL=alignment[10],
            CSTAG=alignment[idx_cstag].replace("cs:Z", ""),
        )
