from __future__ import annotations

import pytest

from csvtag.splitter import split_by_inversion, split_by_tag


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


@pytest.mark.parametrize(
    "csv_tag, expected",
    [
        ("=AA=aa*ga=a=TT", ["=A", "=A", "=a", "=a", "*ga", "=a", "=T", "=T"]),
        ("=A~AA5CC=A", ["=A", "=N", "=N", "=N", "=N", "=N", "=A"]),
        ("=A+TTT=CC-AA=T*AG=TT", ["=A", "+T|+T|+T|=A", "=c", "=a", "-A", "-A", "=T", "*AG", "=T", "=T"]),
    ],
)
def split_by_nucleotide(csv_tag, expected):
    result = list(split_by_inversion(csv_tag))
    assert result == expected, f"Expected {expected}, but got {result}"
