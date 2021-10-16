import os

import pytest

from knotify.rna_analysis import StringAnalyser

GRAMMAR = os.getenv("GRAMMAR", "./libpseudoknot.so")


@pytest.mark.parametrize(
    "sequence, brackets",
    [
        (
            "AUCCUUUUCA",
            [
                "(.....[).]",
                "(...[..).]",
                "(....[)..]",
                "(...[)...]",
                "(....[.).]",
                "(...[.)..]",
            ],
        ),
        (
            "AGAUUGACU",
            [
                "(..[).]..",
                "(.[.)...]",
                "(.[)....]",
            ],
        ),
    ],
)
def test_get_pseudoknots(sequence, brackets):
    s = StringAnalyser(sequence, GRAMMAR, allow_ug=False).get_pseudoknots()

    assert set(a["dot_bracket"] for a in s) == set(brackets)
