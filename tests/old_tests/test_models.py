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
from knotify.models import Pknot


class TestPseudoknot:

    pknot_refer = """...LL.LLL(....RR.RR{..)LLL.LL....}RR.RR....."""
    pknot = {
        "left_core_indices": (22, 9),
        "right_core_indices": (19, 33),
        "left_loop_stems": ([23, 24, 25, 27, 28], [3, 4, 6, 7, 8]),
        "right_loop_stems": ([14, 15, 17, 18], [34, 35, 37, 38]),
    }

    def test_knot_correctly_initialized(self):
        pknot = Pknot(len(self.pknot_refer), **self.pknot)
        assert pknot._repr == self.pknot
        assert pknot.len == len(self.pknot_refer)

    def test_knot_left_outer_part(self):
        knot = Pknot(len(self.pknot_refer), **self.pknot)
        expected = [0, 1, 2, 5]
        assert sorted(knot.get_left_outer_part()) == sorted(expected)

    def test_knot_left_inner_part(self):
        knot = Pknot(len(self.pknot_refer), **self.pknot)
        expected = [10, 11, 12, 13, 16]
        assert sorted(knot.get_left_inner_part()) == sorted(expected)

    def test_knot_right_outer_part(self):
        knot = Pknot(len(self.pknot_refer), **self.pknot)
        expected = [36, 39, 40, 41, 42, 43]
        assert sorted(knot.get_right_outer_part()) == sorted(expected)

    def test_knot_right_inner_part(self):
        knot = Pknot(len(self.pknot_refer), **self.pknot)
        expected = [26, 29, 30, 31, 32]
        assert sorted(knot.get_right_inner_part()) == sorted(expected)

    def test_get_dd_sef(self):
        knot = Pknot(len(self.pknot_refer), **self.pknot)
        expected = [20, 21]
        assert sorted(knot.get_dd_seg()) == sorted(expected)

    def test_get_left_outer_potential(self):
        knot = Pknot(len(self.pknot_refer), **self.pknot)
        expected = (0, 9)
        assert sorted(knot.get_left_outer_potential()) == sorted(expected)

    def test_get_left_inner_potential(self):
        knot = Pknot(len(self.pknot_refer), **self.pknot)
        expected = (23, 33)
        assert knot.get_left_inner_potential() == expected

    def test_get_right_outer_potential(self):
        knot = Pknot(len(self.pknot_refer), **self.pknot)
        expected = (34, len(self.pknot_refer))
        assert knot.get_right_outer_potential() == expected

    def test_right_inner_potential(self):
        knot = Pknot(len(self.pknot_refer), **self.pknot)
        expected = (10, 19)
        assert knot.get_right_inner_potential() == expected
