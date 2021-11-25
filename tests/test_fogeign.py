import pytest

from knotify.foreign import ipknot

IPKNOT_ENERGY = "> case (e=-20.432)\n\n.()."
IPKNOT_E_ENERGY = "> case (e=-1.6e+20)\n\n.().\n"
IPKNOT_NO_ENERGY = "> case\n\n.().\n"


@pytest.mark.parametrize(
    "stdout, result",
    [
        (IPKNOT_ENERGY, {"dot_bracket": ".().", "stems": -1, "energy": -20.432}),
        (IPKNOT_E_ENERGY, {"dot_bracket": ".().", "stems": -1, "energy": -1.6e20}),
        (IPKNOT_NO_ENERGY, {"dot_bracket": ".().", "stems": -1, "energy": 1000}),
    ],
)
def test_ipknot_parse(stdout, result):
    assert ipknot.parse_ipknot_output(stdout) == result
