import pytest
from knotify.pairalign.cpairalign import CPairAlign

from knotify.extensions.skip_final_au import SkipFinalAU

SKIPFINALAU_SO = "./libskipfinalau.so"


@pytest.mark.parametrize(
    "sequence,input,expected",
    [
        (
            "AAGCCUUG",
            ("((([)))]", 2, 0),
            [(".(([)).]", 1, 0)],
        ),
        (
            "AAGACCUUGU",
            ("((([[)))]]", 2, 1),
            [(".(([[)).]]", 1, 1), ("(((.[)))].", 2, 0), (".((.[)).].", 1, 0)],
        ),
        (
            "AAGUCUUA",
            ("((([)))]", 2, 0),
            [(".(([)).]", 1, 0)],
        ),
    ],
)
def test_skip_final_au(sequence: str, input: str, expected: list):
    a = SkipFinalAU(SKIPFINALAU_SO)

    assert set(a.get_candidates(sequence.lower(), *input)) == set(expected)


CPAIRALIGN_SO = "./libcpairalign.so"


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
def test_cpairalign(sequence, core_stems, expected):
    p = CPairAlign(CPAIRALIGN_SO)
    assert p.pairalign(sequence, *core_stems) == [expected]
