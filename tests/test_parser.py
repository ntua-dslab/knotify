#
# Copyright © 2022 Christos Pavlatos, George Rassias, Christos Andrikos,
#                  Evangelos Makris, Aggelos Kolaitis
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the “Software”), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
import os
from typing import Dict, List
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

    result = p.detect_pseudoknots("a{}aaa{}a{}aaa{}aaa".format(A[0], B[0], A[1], B[1]))
    if not combination_has_ug or combination_has_ug and allow_ug:
        assert "1,11,3,1" in result
    if combination_has_ug and not allow_ug:
        assert "1,11,3,1" not in result


@for_each_parser("parser, library_path")
@pytest.mark.parametrize(
    "sequence, result",
    [
        ("caaaaaggaaaaac", ["0,14,5,0"]),
        ("caaaaagagaaaaac", ["0,15,5,1"]),
        ("caaaaagaagaaaaac", ["0,16,5,2"]),
        ("caaaaagaaagaaaaac", ["0,17,5,3"]),
    ],
)
@pytest.mark.parametrize(
    "args",
    [
        {"max_window_size": 100, "min_window_size": 1},
        {
            "max_window_size_ratio": 1.0,
            "min_window_size_ratio": 0,
            "max_window_size": 0,
            "min_window_size": 0,
        },
    ],
)
def test_window_size(
    parser, library_path: str, sequence: str, result: List[str], args: Dict
):
    p = parser(library_path=library_path, max_dd_size=3, allow_ug=False, **args)
    assert p.detect_pseudoknots(sequence) == result
