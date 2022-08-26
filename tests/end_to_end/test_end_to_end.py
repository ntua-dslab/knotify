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

from knotify import knotify

HAIRPIN = os.getenv("HAIRPIN_SO", "./libhairpin.so")


@pytest.mark.parametrize("parser", ["bruteforce", "yaep"])
@pytest.mark.parametrize(
    "name, sequence, expected, config",
    [
        (
            "baseline",
            "GGGAAAUGGACUGAGCGGCGCCGACCGCCAAACAACCGGCA",
            "..............((((.[[[[.))))........]]]].",
            {},
        ),
        (
            "baseline, ug",
            "ACGUGAAGGCUACGAUAGUGCCAG",
            ".((((..[[[)))).....]]]..",
            {"allow_ug": True},
        ),
        (
            "baseline",
            "GGGAAACGGAGUGCGCGGCACCGUCCGCGGAACAAACGGAGAAGGCAGCU",
            ".............(((((..[[[[)))))......]]]]...........",
            {},
        ),
        (
            "baseline",
            "AUCCUUUUCAGUUGGGCCUUCUGGUGAUGUUUCUGGCCACCCAGGAGGUCCUGAGGAAGAGGUGGACGGCCAGAUUGACU",
            ".............(((((((((((.......[[[[[[[..)))))))))))................]]]]]]]......",
            {},
        ),
        (
            "pick smaller dd for same energy, stems",
            "AAAAAACUAAUAGAGGGGGGACUUAGCGCCCCCCAAACCGUAACCCC",
            "..............((((((.....[[[))))))....]]]......",
            {},
        ),
        (
            "dd 7",
            "gggaaacaacaggagggggccacguguggugccguccgcgcccccuauguuguaacagaagcaccacc",
            ".............(((((((.....[[[[[[[.......)))))))..............]]]]]]].",
            {"max_dd_size": 7},
        ),
        (
            "hairpin (disabled)",
            "CGCUCAACUCAUGGAGCGCACGACGGAUCACCAUCGAUUCGACAUGAG",
            "(((((..[[[[[[)))))........................]]]]]]",
            {},
        ),
        (
            "hairpin (disabled)",
            "CGCUCAACUCAUGGAGCGCACGACGGAUCACCAUCGAUUCGACAUGAG",
            "(((((..[[[[[[)))))........................]]]]]]",
            {},
        ),
        (
            "hairpin (enabled)",
            "CGCUCAACUCAUGGAGCGCACGACGGAUCACCAUCGAUUCGACAUGAG",
            "(((((..[[[[[[))))).....((((((......)))))).]]]]]]",
            {"allow_ug": True, "hairpin_grammar": HAIRPIN, "hairpin_allow_ug": True},
        ),
        (
            "hairpin fixes pseudoknot (disabled)",
            "CGCUCAACCUCAGAGCGCAAGAGUCGAACGAAUACAGUUCGACAUGAGG",
            "......................(((((((.....[[))))))).]]...",
            {},
        ),
        (
            "hairpin fixes pseudoknot (enabled)",
            "CGCUCAACCUCAGAGCGCAAGAGUCGAACGAAUACAGUUCGACAUGAGG",
            "(((((..[[[[[))))).....(((((((.......))))))).]]]]]",
            {"hairpin_grammar": HAIRPIN},
        ),
        (
            "hairpin (big)",
            "AGUUCUCCUUAGAGUGUGUGUUGUUUACCCAACGCACUGUCCCUAUGGGGGGCCAACAUAGGUCCA",
            ".(((((((.......((((((((......))))))))....[[[[[)))))))....]]]]]....",
            {"allow_ug": True, "hairpin_grammar": HAIRPIN, "hairpin_allow_ug": True},
        ),
        # TODO(akolaitis): fix this
        # (
        #     "no hairpins",
        #     "UUUAAACUGGUGGGGCAGUGUCUAGGAUUGACGUUAGACACUGCUUUUUGCCCGUUUCAAACAGGUGAAUACAAACCGUCAU",
        #     "...........((((((((((((((...[[[[[)))))))))))))).............................]]]]].",
        #     {"hairpin_grammar": HAIRPIN},
        # ),
        (
            "bulges",
            "GCGUGGAAGCCCUGCCUGGGGUUGAAGCGUUAAAACUUAAUCAGGC",
            "((((...(((((.[[[[[)))))...))))...........]]]]]",
            {"pairalign": "bulges", "max_bulge_size": 3},
        ),
        (
            "bulges fixes pseudoknots (enabled)",
            "GUUUGUUAGUGGCGUGUCCGUCCGCAGCUGGCAAGCGAAUGUAAAGACUGAC",
            "((((((((((.(((.[[[.[[[))).)))))))))).........]]].]]]",
            {"pairalign": "bulges"},
        ),
        (
            "bulges fixes pseudoknots (disabled)",
            "GUUUGUUAGUGGCGUGUCCGUCCGCAGCUGGCAAGCGAAUGUAAAGACUGAC",
            "(((((((((...............[[[)))))))))...........]]]..",
            {},
        ),
    ],
)
def test_knotify(name, parser, sequence, expected, config):
    opts = knotify.new_options()

    opts.parser = parser
    opts.prune_early = True

    for key, value in config.items():
        opts.set_override(key, value)

    algorithm, config = knotify.from_options(opts)
    result = algorithm.get_results(sequence, **config)

    assert result.loc[0].dot_bracket == expected
