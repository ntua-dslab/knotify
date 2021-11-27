import pytest

from knotify.algorithm import ipknot

IPKNOT_ENERGY = "> case (e=-20.432)\n\n.()."
IPKNOT_E_ENERGY = "> case (e=-1.6e+20)\n\n.().\n"
IPKNOT_NO_ENERGY = "> case\n\n.().\n"
IPKNOT_BRACKET_1 = "> case\n\n.([.)].\n"
IPKNOT_BRACKET_2 = "> case\n\n.[(.]).\n"


@pytest.mark.parametrize(
    "stdout, result",
    [
        (IPKNOT_ENERGY, {"dot_bracket": ".().", "stems": 2, "energy": -20.432}),
        (IPKNOT_E_ENERGY, {"dot_bracket": ".().", "stems": 2, "energy": -1.6e20}),
        (IPKNOT_NO_ENERGY, {"dot_bracket": ".().", "stems": 2, "energy": 1000}),
        (IPKNOT_BRACKET_1, {"dot_bracket": ".([.)].", "stems": 4, "energy": 1000}),
        (IPKNOT_BRACKET_2, {"dot_bracket": ".([.)].", "stems": 4, "energy": 1000}),
    ],
)
def test_ipknot_parse(stdout, result):
    assert ipknot.parse_ipknot_output(stdout) == result
