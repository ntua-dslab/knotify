import os

import pytest

from pseudoknot_detector.main import get_results

PSEUDOKNOT = os.getenv("PSEUDOKNOT_SO", "./libpseudoknot.so")
HAIRPIN = os.getenv("HAIRPIN_SO", "./libhairpin.so")


@pytest.mark.parametrize(
    "name, sequence, result, test_params",
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
            "hairpin (enabled)",
            "CGCUCAACUCAUGGAGCGCACGACGGAUCACCAUCGAUUCGACAUGAG",
            "(((((..[[[[[[))))).....((((((......)))))).]]]]]]",
            {"hairpin_grammar": HAIRPIN, "allow_ug": True},
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
            {"hairpin_grammar": HAIRPIN, "allow_ug": True},
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
def test_end_to_end(name, sequence, result, test_params):
    config = {
        "grammar": PSEUDOKNOT,
        "prune_early": True,
    }
    config.update(test_params)

    results = get_results(sequence, **config)
    assert results.loc[0].dot_bracket == result
