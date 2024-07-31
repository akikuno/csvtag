from __future__ import annotations

import pytest

from csvtag.microhomology_handler import _trim_microhomology


@pytest.mark.parametrize("alignments, expected", [
    (
        [{"CSTAG": "=AAATTT"}, {"CSTAG": "=TTTCCC"}],
        [{"CSTAG": "=AAATTT"}, {"CSTAG": "=CCC"}]
    ),
    (
        [{"CSTAG": "=AAAT*gt=T"}, {"CSTAG": "=TTTCCC"}],
        [{"CSTAG": "=AAA"}, {"CSTAG": "=TTTCCC"}]
    ),
    (
        [{"CSTAG": "=GGGAAA"}, {"CSTAG": "=AAACTTT"}, {"CSTAG": "=TTTCCC"}],
        [{"CSTAG": "=GGGAAA"}, {"CSTAG": "=CTTT"}, {"CSTAG": "=CCC"}]
    ),
    (
        [{"CSTAG": "=AAA"}, {"CSTAG": "=AAG"}],
        [{"CSTAG": "=AAA"}, {"CSTAG": "=G"}]
    ),
    (
        [{"CSTAG": "=AAA"}, {"CSTAG": "=TTT"}],
        [{"CSTAG": "=AAA"}, {"CSTAG": "=TTT"}]
    ),
])
def test_trim_microhomology(alignments, expected):
    result = _trim_microhomology(alignments)
    assert result == expected, f"Expected {expected}, but got {result}"
