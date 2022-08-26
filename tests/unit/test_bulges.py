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
    ],
)
def test_bulges(name, sequence, core_stems, config, results):
    cfg = {
        "library_path": BULGES_SO,
        "max_bulge_size": 1,
        "symmetric_bulges": True,
        "min_stems_after_bulge": 1,
    }
    cfg.update(config)

    assert BulgesPairAlign(**cfg).pairalign(sequence, *core_stems) == results
