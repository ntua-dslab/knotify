import pytest

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
