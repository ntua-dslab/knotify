import os
from typing import List
import itertools

import pytest

from tests.utils import for_each_parser


@for_each_parser("parser, library_path")
@pytest.mark.parametrize(
    "dd_size, sequence, result",
    [
        (1, "caaaaaggaaaaac", ["0,14,5,0"]),
        (1, "caaaaagagaaaaac", ["0,15,5,1"]),
        (1, "caaaaagaagaaaaac", [""]),
        (2, "caaaaagaagaaaaac", ["0,16,5,2"]),
        (2, "caaaaagaaagaaaaac", [""]),
        (3, "caaaaagaaagaaaaac", ["0,17,5,3"]),
    ],
)
def test_max_dd_size(
    parser, library_path: str, dd_size: int, sequence: str, result: List[str]
):
    p = parser(library_path=library_path, max_dd_size=dd_size, allow_ug=False)
    assert p.detect_pseudoknots(sequence) == result


@pytest.mark.parametrize(
    "A, B, allow_ug",
    itertools.product(
        ["au", "gc", "ug", "ua", "cg", "gu"],
        ["au", "gc", "ug", "ua", "cg", "gu"],
        [True, False],
    ),
)
@for_each_parser("parser, library_path")
def test_pseudoknot_pairs(parser, library_path, A, B, allow_ug):
    p = parser(library_path=library_path, max_dd_size=1, allow_ug=allow_ug)
    combination_has_ug = any(x in ["gu", "ug"] for x in [A, B])

    result = p.detect_pseudoknots("{}aaa{}a{}aaa{}".format(A[0], B[0], A[1], B[1]))
    if not combination_has_ug or combination_has_ug and allow_ug:
        assert "0,11,3,1" in result
    if combination_has_ug and not allow_ug:
        assert "0,11,3,1" not in result
