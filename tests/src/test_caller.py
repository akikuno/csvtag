from __future__ import annotations

from pathlib import Path

from csvtag.caller import call_csvtag


def test_call_csvtag_one_alignment():
    path_sam = Path("tests/data/one_alignment.sam")
    result = list(call_csvtag(path_sam))
    expected = [{"read1": ["=AAAAA"]}]
    assert result == expected, f"Expected {expected}, but got {result}"


def test_call_csvtag_tow_alignment():
    path_sam = Path("tests/data/two_alignments.sam")
    result = list(call_csvtag(path_sam))
    expected = [{"read1": ["=AAAAA", "=TT*TC=TT"]}]
    assert result == expected, f"Expected {expected}, but got {result}"
