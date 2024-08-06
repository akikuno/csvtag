from __future__ import annotations

import pytest

from csvtag.to_sequence import to_sequence


@pytest.mark.parametrize(
    "csv_tag, expected",
    [
        ("=AA=aa*ga=a=AA", "AAttttAA"),  # 基本的なCSVタグ
        ("=AA=aa+gg=aa=AA", "AAttccttAA"),  # 挿入あり
        ("=AA=aa-gg=aa=AA", "AAttttAA"),  # 欠失あり
        ("=AAGG*CT=AT", "AAGGTAT"),  # 全て大文字のCSVタグ
        ("=tt", "aa"),  # 全て小文字のCSVタグ
        ("=A=a=G", "AtG"),  # 混合した大文字と小文字のCSVタグ
        ("", ""),  # 空のCSVタグ
    ],
)
def test_to_sequence(csv_tag: str, expected: str):
    result = to_sequence(csv_tag)
    assert result == expected, f"Expected {expected}, but got {result}"
