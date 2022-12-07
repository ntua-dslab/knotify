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
from knotify.pairalign.bulges import BulgesPairAlign
from knotify.pairalign.cpairalign import CPairAlign


CPAIRALIGN_SO = "./libcpairalign.so"
BULGES_SO = "./libbulges.so"


@pytest.mark.parametrize(
    "pairalign",
    [
        CPairAlign(CPAIRALIGN_SO),
        BulgesPairAlign(
            max_bulge_size=0,
            min_stems_after_bulge=0,
            symmetric_bulges=True,
            count_stems_from_bulges=False,
            library_path=BULGES_SO,
        ),
    ],
)
@pytest.mark.parametrize(
    "sequence,core_stems,expected",
    [
        ("AAGCCUUG", (2, 6, 0, 0), ("((([)))]", 2, 0)),
        ("AAGACCUUGU", (2, 7, 1, 0), ("((([[)))]]", 2, 1)),
        ("AAGACCUGGU", (2, 7, 1, 0), (".(([[)).]]", 1, 1)),
        ("AAGACCGUGU", (2, 7, 1, 0), ("..([[)..]]", 0, 1)),
        ("AAGACCUGGG", (2, 7, 1, 0), (".((.[)).].", 1, 0)),
        ("AAGACAACUUGG", (2, 9, 1, 2), ("(((.[..)))].", 2, 0)),
        (
            "acgugaaggcuacgauagugccag",
            (4, 16, 4, 0),
            (".((((..[[[)))).....]]]..", 3, 2),
        ),
        (
            "gcguggaagcccugccugggguugaagcguuaaaacuuaaucaggc",
            (11, 31, 5, 0),
            (".......(((((.[[[[[)))))..................]]]]]", 4, 4),
        ),
        (
            "GGGAAACGAGCCAAGUGGCGCCGACCACUUAAAAACACCGGAA",
            (22, 19, 12, 2),
            (".....................((............[..))]..", 1, 0),
        ),
    ],
)
def test_pairalign(pairalign, sequence, core_stems, expected):
    assert pairalign.pairalign(sequence, *core_stems) == [expected]
