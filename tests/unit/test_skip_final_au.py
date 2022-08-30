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
import pytest

from knotify.extensions.skip_final_au import SkipFinalAU

SKIPFINALAU_SO = "./libskipfinalau.so"


@pytest.mark.parametrize(
    "sequence,input,expected",
    [
        (
            "AAGCCUUG",
            ("((([)))]", 2, 0),
            [(".(([)).]", 1, 0)],
        ),
        (
            "AAGACCUUGU",
            ("((([[)))]]", 2, 1),
            [(".(([[)).]]", 1, 1), ("(((.[)))].", 2, 0), (".((.[)).].", 1, 0)],
        ),
        (
            "AAGUCUUA",
            ("((([)))]", 2, 0),
            [(".(([)).]", 1, 0)],
        ),
    ],
)
def test_skip_final_au(sequence: str, input: str, expected: list):
    a = SkipFinalAU(SKIPFINALAU_SO)

    assert set(a.get_candidates(sequence.lower(), *input)) == set(expected)
