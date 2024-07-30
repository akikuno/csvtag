from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from csvtag.sam_handler import extract_alignment, read_sam, remove_resequence

# 入力: samのalignment: Iterator[str]
# 出力: {QNAME: csvtag}: Iterator[dict[str, str]]
# Inversionはリファレンスに合わせる（＝revcompする）


def call_csvtag(path_sam: str | Path) -> dict[str, str]:
    alignments: Iterator[dict[str, str | int]] = extract_alignment(read_sam(path_sam))
    alignments = remove_resequence(alignments)
