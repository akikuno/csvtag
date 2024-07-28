from __future__ import annotations

from src.csvtag.split import split
import pytest

@pytest.mark.parametrize("csv_tag, expected", [
    (":4*AG:3", [':4', '*AG', ':3']),
    ("=AA=aa=TT", ['=AA', '=aa', '=TT']),
    ("=AA=aa*ga=a=TT", ['=AA', '=aa', '*ga', '=a', '=TT']),
    (":10+AGC:3", [':10', '+AGC', ':3']),
    ("~ACGT12ACGT", ['~ACGT12ACGT']),
    ("=NN:5-TC", ['=NN', ':5', '-TC']),
])
def test_split(csv_tag, expected):
    result = split(csv_tag)
    assert result == expected, f"Expected {expected}, but got {result}"
