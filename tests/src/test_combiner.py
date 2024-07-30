from __future__ import annotations

import pytest

from csvtag.combiner import combine_splitted_csv_tag


@pytest.mark.parametrize(
    "splitted_csv_tag, expected",
    [
        (iter(["=A", "+T|+T|+T|=C", "=C", "-A", "-A", "=T", "*AG", "=T", "=T"]), "=A+TTT=CC-AA=T*AG=TT"),
        (iter(["=A", "=A", "=A", "=N", "=N", "=N", "=N", "=N", "=C", "=C", "=A"]), "=AAANNNNNCCA"),
        (iter(["*AG", "*T", "*C", "*G"]), "*AG*T*C*G"),
        (iter(["=G", "-A", "+T", "=C"]), "=G-A+T=C"),
    ],
)
def test_combine_splitted_csv_tag(splitted_csv_tag, expected):
    result = combine_splitted_csv_tag(splitted_csv_tag)
    assert result == expected, f"Expected {expected}, but got {result}"
