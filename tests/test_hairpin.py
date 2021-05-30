import os

import pytest
import pandas as pd

from pseudoknot_detector.hairpin import (
    HairpinDetector,
    get_loop_indices,
    find_hairpins,
    dot_bracket_to_record,
)


@pytest.mark.parametrize(
    "bracket,indices",
    [
        ("([)]", (1, 1, 3, 3)),
        ("([)......]", (1, 1, 3, 9)),
        ("(......[)......]", (1, 7, 9, 15)),
        ("(......[.).....]", (1, 7, 10, 15)),
        (".(((((((................[[[[))))))).................]]]].", (8, 24, 35, 52)),
    ],
)
def test_get_loop_indices(bracket: str, indices: tuple):
    lstart, lend, rstart, rend = get_loop_indices(bracket)
    assert (lstart, lend, rstart, rend) == indices

    assert set(bracket[lstart:lend]) <= {"."}
    assert set(bracket[rstart:rend]) <= {"."}


HAIRPIN = os.getenv("HAIRPIN_SO", "./libhairpin.so")


@pytest.mark.parametrize(
    "sequence, result",
    [
        ("UAAAUCGUCGUGAGACG", ".....((((....))))"),
        ("AAGAUGGUUUAAACCA", "....((((....))))"),
        ("AUAAACCUCAAUAAUCCCUCGUUAGGAGGGA", "..............((((((.....))))))"),
        ("CAAGAGUCGAACGAAUACAGUUCGACA", ".....(((((((.......)))))))."),
        ("UUAAAUAGUUUAUCUGUUAAGAUAAAC", ".......((((((((....))))))))"),
    ],
)
def test_hairpin_detector(sequence: str, result: str):
    detector = HairpinDetector(grammar=HAIRPIN, allow_ug=False)
    assert result in detector.detect_hairpins(sequence)


@pytest.mark.parametrize(
    "sequence, result, result_ug",
    [
        (
            "UUAGAGUGUGUGUUGUUUACCCAACGCACUGUC",
            "....((((..................))))...",
            ".......((((((((......))))))))....",
        ),
        (
            "CACGACGGAUCACCAUCGAUUCGA",
            "..(((...............))).",
            ".....((((((......)))))).",
        ),
    ],
)
def test_hairpin_detector_allow_ug(sequence: str, result: str, result_ug: str):
    res = HairpinDetector(grammar=HAIRPIN, allow_ug=False).detect_hairpins(sequence)
    res_ug = HairpinDetector(grammar=HAIRPIN, allow_ug=True).detect_hairpins(sequence)

    assert result in res
    assert result in res_ug
    assert result_ug not in res
    assert result_ug in res_ug


@pytest.mark.parametrize(
    "dot_bracket, left_loop_stems, right_loop_stems, dd",
    [
        ("([)]", 0, 0, 0),
        ("((((([)))))]", 4, 0, 0),
        ("((((([..)))))]", 4, 0, 2),
        (".((((([..)))))]", 4, 0, 2),
        (".((((([[[[..)))))]]]]", 4, 3, 2),
        (".((((([[[[..)))))]]]].......", 4, 3, 2),
        (".((((([[[[..)))))..........]]]].......", 4, 3, 2),
        (".(((((..............[[[[..)))))..........]]]].......", 4, 3, 2),
    ],
)
def test_dot_bracket_to_record(dot_bracket, left_loop_stems, right_loop_stems, dd):
    assert dot_bracket_to_record(dot_bracket) == (left_loop_stems, right_loop_stems, dd)


def test_find_hairpins():
    sequence = "CGGUAGAAAAGAUGGUUUAAACCACGCCUUCUACCAAGUUAGUAAAUAAAUAGGCGG"
    data = pd.DataFrame.from_records(
        [
            {
                "dot_bracket": ".(((((((................[[[[))))))).................]]]].",
            },
            {
                "dot_bracket": ".((((((.................[[[[[))))))................]]]]].",
            },
        ]
    )

    data_with_hairpins = find_hairpins(sequence, data, HAIRPIN, min_stems=4, min_size=3)

    assert set(data_with_hairpins["dot_bracket"].unique()) == set(
        [
            ".(((((((................[[[[))))))).................]]]].",
            ".(((((((....((((....))))[[[[))))))).................]]]].",
            ".((((((.................[[[[[))))))................]]]]].",
            ".((((((.....((((....))))[[[[[))))))................]]]]].",
        ]
    )


@pytest.mark.parametrize(
    "name, sequence, config, results",
    [
        (
            "normal",
            "augc",
            {"max_per_loop": 1},
            ["()..", "..()"],
        ),
        (
            "two hairpins",
            "augc",
            {"max_per_loop": 2},
            ["()..", "..()", "()()"],
        ),
        (
            "two hairpins",
            "augc",
            {"max_per_loop": 2, "min_stems": 2},
            ["()()"],
        ),
        (
            "bulge",
            "aaugu",
            {"max_bulge": 1, "min_stems": 2},
            ["(().)"],
        ),
        (
            "big bulge",
            "aaugggu",
            {"max_bulge": 3, "min_stems": 2},
            ["(()...)"],
        ),
        (
            "allow ug",
            "gu",
            {"allow_ug": True},
            ["()"],
        ),
    ],
)
def test_parameters(name: str, sequence: str, config: dict, results: list):
    test_config = {"allow_ug": False, "min_size": 0, "min_stems": 1}
    test_config.update(config)
    detector = HairpinDetector(HAIRPIN, **test_config)
    assert set(detector.detect_hairpins(sequence)) == set(results)
