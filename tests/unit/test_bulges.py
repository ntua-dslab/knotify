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

BULGES_SO = "./libbulges.so"


@pytest.mark.parametrize(
    "name, sequence, core_stems, config, results",
    [
        (
            "simple",
            "GCGUGGAAGCCCUGCCUGGGGUUGAAGCGUUAAAACUUAAUCAGGC",
            (11, 31, 5, 0),
            {"max_bulge_size": 3, "min_stems_after_bulge": 4, "symmetric_bulges": True},
            [
                (".......(((((.[[[[[)))))..................]]]]]", 4, 4),
                ("((((...(((((.[[[[[)))))...))))...........]]]]]", 8, 4),
            ],
        ),
        (
            "maxsize",
            "GCGUGGAAGCCCUGCCUGGGGUUGAAGCGUUAAAACUUAAUCAGGC",
            (11, 31, 5, 0),
            {"max_bulge_size": 2, "min_stems_after_bulge": 4, "symmetric_bulges": True},
            [
                (".......(((((.[[[[[)))))..................]]]]]", 4, 4),
            ],
        ),
        (
            "minstems",
            "GCGUGGAAGCCCUGCCUGGGGUUGAAGCGUUAAAACUUAAUCAGGC",
            (11, 31, 5, 0),
            {"max_bulge_size": 3, "min_stems_after_bulge": 5, "symmetric_bulges": True},
            [
                (".......(((((.[[[[[)))))..................]]]]]", 4, 4),
            ],
        ),
        (
            "assymetric",
            "GCGUGGAAAGCCCUGCCUGGGGUUGAAGCGUUAAAACUUAAUCAGGC",
            (12, 31, 5, 0),
            {
                "max_bulge_size": 4,
                "min_stems_after_bulge": 4,
                "symmetric_bulges": False,
            },
            [
                ("........(((((.[[[[[)))))..................]]]]]", 4, 4),
                ("((((....(((((.[[[[[)))))...))))...........]]]]]", 8, 4),
            ],
        ),
        (
            "from-dataset",
            "GUGGCAGUAAGCCUGGGAAUGGGGGCGACCCAGGCGUAUGAACAUAGUGUAACGCUCCCC",
            (16, 37, 9, 1),
            {"max_bulge_size": 1, "min_stems_after_bulge": 3, "symmetric_bulges": True},
            [
                ("..........(((((((...[[[[[[[.))))))).................]]]]]]].", 6, 6),
                ("......(((.(((((((...[[[[[[[.))))))).))).............]]]]]]].", 9, 6),
            ],
        ),
        (
            "right",
            "AGGGGGGGUUAAUUUUUU",
            (0, 11, 7, 0),
            {"max_bulge_size": 1, "symmetric_bulges": False},
            [
                ("(.......[)].......", 0, 0),
                ("([[[[[[.[)].]]]]]]", 0, 6),
            ],
        ),
        (
            "both-sides",
            "CAAAACAGGGGGGGUUCUUUUAAUUUUUU",
            (6, 16, 7, 0),
            {"max_bulge_size": 1, "symmetric_bulges": False},
            [
                ("......(.......[).....].......", 0, 0),
                (".((((.(.......[).))))].......", 4, 0),
                ("......([[[[[[.[).....].]]]]]]", 0, 6),
                (".((((.([[[[[[.[).))))].]]]]]]", 4, 6),
            ],
        ),
        (
            "both-sides-no-count",
            "CAAAACAGGGGGGGUUCUUUUAAUUUUUU",
            (6, 16, 7, 0),
            {
                "max_bulge_size": 1,
                "symmetric_bulges": False,
                "count_stems_from_bulges": False,
            },
            [
                ("......(.......[).....].......", 0, 0),
                (".((((.(.......[).))))].......", 0, 0),
                ("......([[[[[[.[).....].]]]]]]", 0, 0),
                (".((((.([[[[[[.[).))))].]]]]]]", 0, 0),
            ],
        ),
    ],
)
def test_bulges(name, sequence, core_stems, config, results):
    cfg = {
        "library_path": BULGES_SO,
        "max_bulge_size": 1,
        "symmetric_bulges": True,
        "min_stems_after_bulge": 1,
        "count_stems_from_bulges": True,
    }
    cfg.update(config)

    assert BulgesPairAlign(**cfg).pairalign(sequence, *core_stems) == results
