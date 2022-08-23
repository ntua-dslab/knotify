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

import pytest

from knotify.algorithm.knotify import Knotify
from knotify.pairalign.cpairalign import CPairAlign
from knotify.pairalign.skip_final_au import SkipFinalAU
from tests.test_pairalign import CPAIRALIGN_SO, SKIPFINALAU_SO
from tests.utils import for_each_parser

HAIRPIN = os.getenv("HAIRPIN_SO", "./libhairpin.so")


@pytest.mark.parametrize(
    "name, sequence, result, parser_params, test_params",
    [
        (
            "baseline",
            "GGGAAAUGGACUGAGCGGCGCCGACCGCCAAACAACCGGCA",
            "..............((((.[[[[.))))........]]]].",
            {},
            {},
        ),
        (
            "baseline, ug",
            "ACGUGAAGGCUACGAUAGUGCCAG",
            ".((((..[[[)))).....]]]..",
            {"allow_ug": True},
            {},
        ),
        (
            "baseline",
            "GGGAAACGGAGUGCGCGGCACCGUCCGCGGAACAAACGGAGAAGGCAGCU",
            ".............(((((..[[[[)))))......]]]]...........",
            {},
            {},
        ),
        (
            "baseline",
            "AUCCUUUUCAGUUGGGCCUUCUGGUGAUGUUUCUGGCCACCCAGGAGGUCCUGAGGAAGAGGUGGACGGCCAGAUUGACU",
            ".............(((((((((((.......[[[[[[[..)))))))))))................]]]]]]]......",
            {},
            {},
        ),
        (
            "pick smaller dd for same energy, stems",
            "AAAAAACUAAUAGAGGGGGGACUUAGCGCCCCCCAAACCGUAACCCC",
            "..............((((((.....[[[))))))....]]]......",
            {},
            {},
        ),
        (
            "dd 7",
            "gggaaacaacaggagggggccacguguggugccguccgcgcccccuauguuguaacagaagcaccacc",
            ".............(((((((.....[[[[[[[.......)))))))..............]]]]]]].",
            {"max_dd_size": 7},
            {},
        ),
        (
            "hairpin (disabled)",
            "CGCUCAACUCAUGGAGCGCACGACGGAUCACCAUCGAUUCGACAUGAG",
            "(((((..[[[[[[)))))........................]]]]]]",
            {},
            {},
        ),
        (
            "hairpin (enabled)",
            "CGCUCAACUCAUGGAGCGCACGACGGAUCACCAUCGAUUCGACAUGAG",
            "(((((..[[[[[[))))).....((((((......)))))).]]]]]]",
            {"allow_ug": True},
            {"hairpin_grammar": HAIRPIN, "hairpin_allow_ug": True},
        ),
        (
            "hairpin fixes pseudoknot (disabled)",
            "CGCUCAACCUCAGAGCGCAAGAGUCGAACGAAUACAGUUCGACAUGAGG",
            "......................(((((((.....[[))))))).]]...",
            {},
            {},
        ),
        (
            "hairpin fixes pseudoknot (enabled)",
            "CGCUCAACCUCAGAGCGCAAGAGUCGAACGAAUACAGUUCGACAUGAGG",
            "(((((..[[[[[))))).....(((((((.......))))))).]]]]]",
            {},
            {"hairpin_grammar": HAIRPIN},
        ),
        (
            "hairpin (big)",
            "AGUUCUCCUUAGAGUGUGUGUUGUUUACCCAACGCACUGUCCCUAUGGGGGGCCAACAUAGGUCCA",
            ".(((((((.......((((((((......))))))))....[[[[[)))))))....]]]]]....",
            {"allow_ug": True},
            {"hairpin_grammar": HAIRPIN, "hairpin_allow_ug": True},
        ),
        # TODO(akolaitis): fix this
        # (
        #     "no hairpins",
        #     "UUUAAACUGGUGGGGCAGUGUCUAGGAUUGACGUUAGACACUGCUUUUUGCCCGUUUCAAACAGGUGAAUACAAACCGUCAU",
        #     "...........((((((((((((((...[[[[[)))))))))))))).............................]]]]].",
        #     {"hairpin_grammar": HAIRPIN},
        # ),
    ],
)
@for_each_parser("parser, library_path")
def test_end_to_end(
    parser, library_path, name, sequence, result, parser_params, test_params
):
    config = {
        "parser": parser(library_path, **parser_params),
        "pairalign": CPairAlign(CPAIRALIGN_SO).pairalign,
        "prune_early": True,
    }
    config.update(test_params)

    results = Knotify().get_results(sequence, **config)
    assert results.loc[0].dot_bracket == result


@pytest.mark.parametrize(
    "name, sequence, candidate, test_params",
    [
        (
            "without final au",
            "GGGAAACGAGCCAAGUGGCGCCGACCACUUAAAAACACCGGAA",
            ".............(((((..[[[.))))).........]]]..",
            {"skip_final_au": SkipFinalAU(SKIPFINALAU_SO).pairalign},
        ),
        *[
            (
                "without final au (both sides)",
                "GGGAAACGAGCCAAGUGGCUCCGACCACUUAAAAACACCGGAA",
                candidate,
                {
                    "skip_final_au": SkipFinalAU(SKIPFINALAU_SO).pairalign,
                    "max_stem_allow_smaller": 3,
                },
            )
            for candidate in [
                ".............(((((..[[[.))))).........]]]..",  # drop both
                ".............(((((.[[[[.))))).........]]]].",  # drop left
                "............((((((..[[[.))))))........]]]..",  # drop right
                "............((((((.[[[[.))))))........]]]].",  # include all
            ]
        ],
    ],
)
@for_each_parser("parser, library_path")
def test_exploration(parser, library_path, name, sequence, candidate, test_params):
    config = {
        "parser": parser(library_path),
        "pairalign": CPairAlign(CPAIRALIGN_SO).pairalign,
        "prune_early": True,
    }
    config.update(test_params)

    assert candidate in Knotify().get_results(sequence, **config).dot_bracket.unique()
