from __future__ import annotations

import pytest

from csvtag.revcomp import revcomp


@pytest.mark.parametrize(
    "input_str, expected",
    [
        ("=AA=aa*ga=a=AA", "=TT=t*ct=tt=TT"),
        ("=AA+accc=CC", "=GG+gggt=TT"),
        ("=AA~AC10TG=CC", "=GG~CA10GT=TT"),
        ("", ""),
        ("=N", "=N"),
    ],
)
def test_revcomp(input_str, expected):
    result = revcomp(input_str)
    assert result == expected, f"Expected {expected}, but got {result}"
