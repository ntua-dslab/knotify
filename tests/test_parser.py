import os
from typing import List
import itertools

import pytest

from knotify.rna_parser import PseudoknotDetector

PSEUDOKNOT = os.getenv("PSEUDOKNOT_SO", "./libpseudoknot.so")


@pytest.mark.parametrize(
    "dd_size, sequence, result",
    [
        (1, "caaaaaggaaaaac", ["5, 0"]),
        (1, "caaaaagagaaaaac", ["5, 1"]),
        (1, "caaaaagaagaaaaac", [""]),
        (2, "caaaaagaagaaaaac", ["5, 2"]),
        (2, "caaaaagaaagaaaaac", [""]),
        (3, "caaaaagaaagaaaaac", ["5, 3"]),
    ],
)
def test_max_dd_size(dd_size: int, sequence: str, result: List[str]):
    parser = PseudoknotDetector(grammar=PSEUDOKNOT, max_dd_size=dd_size, allow_ug=False)
    assert parser.detect_pseudoknots(sequence) == result


@pytest.mark.parametrize(
    "A, B, allow_ug",
    itertools.product(
        ["au", "gc", "ug", "ua", "cg", "gu"],
        ["au", "gc", "ug", "ua", "cg", "gu"],
        [True, False],
    ),
)
def test_pseudoknot_pairs(A, B, allow_ug):
    parser = PseudoknotDetector(grammar=PSEUDOKNOT, max_dd_size=4, allow_ug=allow_ug)
    combination_has_ug = any(x in ["gu", "ug"] for x in [A, B])

    result = parser.detect_pseudoknots(
        "".join(
            [
                A[0],  # left core outer
                "augcaugc",  # left loop
                B[0],  # right core inner
                "augc",  # dd
                A[1],  # left core inner
                "augcaugc",  # right loop
                B[1],  # right core outer
            ]
        )
    )
    if not combination_has_ug or combination_has_ug and allow_ug:
        assert "8, 4" in result
    if combination_has_ug and not allow_ug:
        assert "8, 4" not in result
