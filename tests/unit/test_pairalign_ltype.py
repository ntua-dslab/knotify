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
from knotify.pairalign.cpairalign import CLTypePairAlign


CPAIRALIGN_LTYPE_SO = "./libcpairalign_ltype.so"


@pytest.mark.parametrize(
    "sequence,core_stems,expected",
    [
        ("caaaacaaggaaagaaac", (0, 18, 4, 2, 3, 0), ("(....{..[)...}...]", 0, 0, 0)),
        ("AUAGUAUGUAUAUAGU", (2, 12, 1, 1, 2, 1), ("((({{.[.)))}}]..", 2, 1, 0)),
        ("AUAGUCUGUAUAUAGU", (2, 12, 1, 1, 2, 1), ("((({{[[.)))}}]].", 2, 1, 1)),
        ("AUAGUCUGUAAAUAGU", (2, 12, 1, 1, 2, 1), (".(({{[[.)).}}]].", 1, 1, 1)),
    ],
)
def test_pairalign(sequence, core_stems, expected):
    pairalign = CLTypePairAlign(library_path=CPAIRALIGN_LTYPE_SO)

    assert pairalign.pairalign(sequence, *core_stems) == [expected]
