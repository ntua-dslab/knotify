import os

import pytest

from knotify.main import get_results
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
@for_each_parser("parser, library_path")
def test_exploration(parser, library_path, name, sequence, candidate, test_params):
    config = {
        "parser": parser(library_path),
        "prune_early": True,
    }
    config.update(test_params)

    assert candidate in get_results(sequence, **config).dot_bracket.unique()
