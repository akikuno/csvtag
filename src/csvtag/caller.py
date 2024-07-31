from __future__ import annotations

from collections.abc import Iterator
from itertools import groupby
from pathlib import Path

from csvtag.microhomology_handler import remove_microhomology
from csvtag.sam_handler import extract_alignment, read_sam, remove_overlapped_alignments

# 入力: samのalignment: Iterator[str]
# 出力: {QNAME: csvtag}: Iterator[dict[str, str]]
# Inversionはリファレンスに合わせる（＝revcompする）


def call_csvtag(path_sam: str | Path) -> dict[str, str]:
    alignments: Iterator[dict[str, str | int]] = extract_alignment(read_sam(path_sam))
    # バグを未然に防ぐため、マイクロホモロジー以外の重複を除く
    alignments = remove_overlapped_alignments(alignments)
    # cstagにあるマイクロホモロジーをトリムする
    alignments = remove_microhomology(alignments)
    # 逆位を検出し、revcompする
    alignments = list(alignments)
    alignments.sort(key=lambda x: (x["QNAME"], x["POS"]))
    for _, alignments_grouped in groupby(alignments, key=lambda x: x["QNAME"]):
        alignments_grouped = list(alignments_grouped)
        if len(alignments_grouped) == 1:
            yield from alignments_grouped
            continue
        # 逆位の検出：
    # {QNAME: CSVTAG}の辞書を作成する
