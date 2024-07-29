from __future__ import annotations

import pytest

from csvtag.spliter import split_by_inversion, split_by_tag


@pytest.mark.parametrize(
    "csv_tag, expected",
    [
        (":4*AG:3", [":4", "*AG", ":3"]),
        ("=AA=aa=TT", ["=AA", "=aa", "=TT"]),
        ("=AA=aa*ga=a=TT", ["=AA", "=aa", "*ga", "=a", "=TT"]),
        (":10+AGC:3", [":10", "+AGC", ":3"]),
        ("~ACGT12ACGT", ["~ACGT12ACGT"]),
        ("=NN:5-TC", ["=NN", ":5", "-TC"]),
    ],
)
def test_split_by_tag(csv_tag, expected):
    result = list(split_by_tag(csv_tag))
    assert result == expected, f"Expected {expected}, but got {result}"


@pytest.mark.parametrize(
    "csv_tag, expected",
    [
        ("=AA=aa*ga=a=TT", ["=AA", "=aa*ga=a", "=TT"]),
        ("=AA=aa=TT", ["=AA", "=aa", "=TT"]),
        ("=AA=tt*cc=GG", ["=AA", "=tt*cc", "=GG"]),
        ("=NN=nn=MM", ["=NN", "=nn", "=MM"]),
        ("=CC=cc*gg=AA", ["=CC", "=cc*gg", "=AA"]),
    ],
)
def test_split_by_case(csv_tag, expected):
    result = list(split_by_inversion(csv_tag))
    assert result == expected, f"Expected {expected}, but got {result}"
