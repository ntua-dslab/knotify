import os

import pytest

from knotify.energy.pkenergy import PKEnergy

PKENERGY_SO = os.getenv("PKENERGY_SO", "./libpkenergy.so")
PKENERGY_PARAMS = os.getenv("PKENERGY_PARAMS", ".pkenergy/hotknots/params")


@pytest.mark.parametrize(
    "sequence, dot_bracket, result",
    [
        (
            "AAUGCAACUUUUAAAUAGUUUAUCUGUUAAGAUAAACCACCUAGGUUGCAUAUAUAAAAAAUAAAAGGUGCC",
            ".(((((((((...........................[[[[[)))))))))..............]]]]]..",
            -6.985999584197998,
        ),
    ],
)
def test_pkenergy(sequence: str, dot_bracket: str, result: float):
    e = PKEnergy(PKENERGY_SO, PKENERGY_PARAMS)
    assert e.energy_eval(sequence, dot_bracket) == result
