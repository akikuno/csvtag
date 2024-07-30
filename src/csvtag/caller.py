from __future__ import annotations

# 入力: samのalignment: Iterator[str]
# 出力: {QNAME: csvtag}: Iterator[dict[str, str]]
# Inversionはリファレンスに合わせる（＝revcompする）
