from __future__ import annotations

import re
from collections.abc import Iterator
from pathlib import Path

###########################################################
# Utility functions
###########################################################


def read_sam(path_of_sam: str | Path) -> Iterator[list[str]]:
    with open(path_of_sam) as f:
        for line in f:
            yield line.strip().split("\t")


def is_forward_strand(flag: int) -> bool:
    return (flag & 0x10) == 0


def split_cigar(cigar: str) -> Iterator[str]:
    cigar_iter = iter(re.split(r"([MIDNSHPX=])", cigar))
    return (i + op for i, op in zip(cigar_iter, cigar_iter))


def calculate_alignment_length(cigar: str) -> int:
    return sum(int(c[:-1]) for c in split_cigar(cigar) if c[-1] in "MDN=X")


def trim_softclip(qual: str, cigar: str) -> str:

    cigar_split = list(split_cigar(cigar))

    if cigar_split[0].endswith("S"):
        softclip_length = int(cigar_split[0][:-1])
        qual = qual[softclip_length:]

    if cigar_split[-1].endswith("S"):
        softclip_length = int(cigar_split[-1][:-1])
        qual = qual[:-softclip_length]

    return qual


###########################################################
# Format headers and alignments
###########################################################


def extract_sqheaders(sam: Iterator[list[str]]) -> dict[str, int]:
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
            CSTAG=alignment[idx_cstag].replace("cs:Z:", ""),
        )
