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
