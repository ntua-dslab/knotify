import os

import pytest

from knotify.main import get_results

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


@pytest.mark.parametrize(
    "name, sequence, candidate, test_params",
    [
        (
            "without final au",
            "GGGAAACGAGCCAAGUGGCGCCGACCACUUAAAAACACCGGAA",
            ".............(((((..[[[.))))).........]]]..",
            {"allow_skip_final_au": True},
        ),
        *[
            (
                "without final au (both sides)",
                "GGGAAACGAGCCAAGUGGCUCCGACCACUUAAAAACACCGGAA",
                candidate,
                {"allow_skip_final_au": True, "max_stem_allow_smaller": 3},
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
def test_exploration(name, sequence, candidate, test_params):
    config = {
        "grammar": PSEUDOKNOT,
        "prune_early": True,
    }
    config.update(test_params)

    assert candidate in get_results(sequence, **config).dot_bracket.unique()
