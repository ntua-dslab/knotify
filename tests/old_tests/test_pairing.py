#
# Copyright © 2021 Christos Pavlatos, George Rassias, Christos Andrikos,
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
import pytest

from knotify.models import Pknot, INNER, OUTER

from knotify.pairing import (
    #                             pairalign,
    #                             choose,
    get_left_stem_aligned_indices,
    get_right_stem_aligned_indices,
    #                             stem_align,
    #                             custom_pairalign
)

from knotify.pairalign.custom_align import custom_pairalign

pairalign = custom_pairalign


@pytest.mark.parametrize(
    "A, B, idxA, idxB, matchA, matchB",
    [
        ("acg", "cgu", 0, 2, "acg", "cgu"),
        ("", "cgu", 0, -1, "", ""),
        ("acg", "", 3, -1, "", ""),
        ("", "", 0, -1, "", ""),
    ],
)
def test_pairalign(A: str, B: str, idxA: int, idxB: int, matchA: str, matchB: str):
    assert pairalign(A, B) == ((idxA, idxB), (matchA, matchB))


def test_get_left_stem_aligned_indices():
    left, right = "acgu", "acgauagu"
    string = left + right

    free_idx, bound_idx = get_left_stem_aligned_indices(string, (4, 12), (0, 4))

    assert "".join([string[i] for i in sorted(bound_idx[0])]) == "acg"
    assert "".join([string[i] for i in sorted(bound_idx[1])]) == "cgu"
    assert "".join(string[i] for i in sorted(free_idx[0])) == "uagu"
    assert "".join(string[i] for i in sorted(free_idx[1])) == "a"


def test_get_right_stem_aligned_indices():
    left, right = "acgauagu", "acgu"
    string = left + right

    free_idx, bound_idx = get_right_stem_aligned_indices(string, (8, 12), (0, 8))

    assert "".join([string[i] for i in sorted(bound_idx[0])]) == "cgu"
    assert "".join([string[i] for i in sorted(bound_idx[1])]) == "acg"
    assert "".join(string[i] for i in sorted(free_idx[0])) == "a"
    assert "".join(string[i] for i in sorted(free_idx[1])) == "uagu"


class TestLeftRightStems:
    def test_get_left_stem_aligned_indices(self):
        input_string = "acgugaaggcuacgauagugccag"
        ground_truth = ".LLL(..RR[)LLL.....]RR.."
        expected_outer_part = [1, 2, 3]
        expected_inner_part = [11, 12, 13]
        left_core_indices = (ground_truth.index(")"), ground_truth.index("("))
        right_core_indices = (ground_truth.index("["), ground_truth.index("]"))

        outer_potential_indices = (0, left_core_indices[OUTER])
        inner_potential_indices = (
            left_core_indices[INNER] + 1,
            right_core_indices[OUTER],
        )
        _, left_stem_indices = get_left_stem_aligned_indices(
            input_string, inner_potential_indices, outer_potential_indices
        )

        assert sorted(left_stem_indices[OUTER]) == sorted(expected_outer_part)
        assert sorted(left_stem_indices[INNER]) == sorted(expected_inner_part)

    def test_get_right_stem_alignment_indices(self):
        input_string = "acgugaaggcuacgauagugccag"
        ground_truth = ".LLL(..RR[)LLL.....]RR.."
        expected_outer_part = [20, 21]
        expected_inner_part = [7, 8]
        left_core_indices = (ground_truth.index(")"), ground_truth.index("("))
        right_core_indices = (ground_truth.index("["), ground_truth.index("]"))

        outer_potential_indices = (
            right_core_indices[OUTER] + 1,
            len(input_string),
        )
        inner_potential_indices = (
            left_core_indices[OUTER] + 1,
            right_core_indices[INNER],
        )
        _, right_stem_indices = get_right_stem_aligned_indices(
            input_string, inner_potential_indices, outer_potential_indices
        )
        assert sorted(right_stem_indices[OUTER]) == sorted(expected_outer_part)
        assert sorted(right_stem_indices[INNER]) == sorted(expected_inner_part)

    def test_get_right_stem_alignment_indices_7(self):
        input_string = "gcguggaagcccugccugggguugaagcguuaaaacuuaaucaggc"
        ground_truth = ".......LLLL(.RRRR{)LLLL..................}RRRR"
        expected_outer_part = [42, 43, 44, 45]
        expected_inner_part = [13, 14, 15, 16]
        left_core_indices = (ground_truth.index(")"), ground_truth.index("("))
        right_core_indices = (ground_truth.index("{"), ground_truth.index("}"))

        outer_potential_indices = (
            right_core_indices[OUTER] + 1,
            len(input_string),
        )
        inner_potential_indices = (
            left_core_indices[OUTER] + 1,
            right_core_indices[INNER],
        )
        _, right_stem_indices = get_right_stem_aligned_indices(
            input_string, inner_potential_indices, outer_potential_indices
        )
        assert sorted(right_stem_indices[OUTER]) == sorted(expected_outer_part)
        assert sorted(right_stem_indices[INNER]) == sorted(expected_inner_part)

    def test_get_left_stem_aligned_indices_7(self):
        input_string = "gcguggaagcccugccugggguugaagcguuaaaacuuaaucaggc"
        ground_truth = ".......LLLL(.RRRR[)LLLL..................]RRRR"
        expected_outer_part = [7, 8, 9, 10]
        expected_inner_part = [19, 20, 21, 22]
        left_core_indices = (ground_truth.index(")"), ground_truth.index("("))
        right_core_indices = (ground_truth.index("["), ground_truth.index("]"))

        outer_potential_indices = (0, left_core_indices[OUTER])
        inner_potential_indices = (
            left_core_indices[INNER] + 1,
            right_core_indices[OUTER],
        )
        _, left_stem_indices = get_left_stem_aligned_indices(
            input_string, inner_potential_indices, outer_potential_indices
        )

        assert sorted(left_stem_indices[OUTER]) == sorted(expected_outer_part)
        assert sorted(left_stem_indices[INNER]) == sorted(expected_inner_part)
