import pytest

from knotify.algorithm import ipknot
from knotify.algorithm import knotty

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


KNOTTY_OUTPUT = "Seq: AUGAAAG\nRES: .([.)].  -3.2\n"
KNOTTY_OUTPUT_ALL = "Seq: AUGAAAAG\nRES: .([{})].  -3.21\n"


@pytest.mark.parametrize(
    "stdout, result",
    [
        (KNOTTY_OUTPUT, {"dot_bracket": ".([.)].", "stems": 4, "energy": -3.2}),
        (KNOTTY_OUTPUT_ALL, {"dot_bracket": ".([{})].", "stems": 6, "energy": -3.21}),
    ],
)
def test_knotty_parse(stdout, result):
    assert knotty.parse_knotty_output(stdout) == result


@pytest.mark.parametrize(
    "algorithm, config",
    [
        (knotty.Knotty(), {"knotty_executable": "./.knotty/knotty"}),
        (ipknot.IPKnot(), {"ipknot_executable": "./.ipknot/ipknot"}),
    ],
)
def test_foreign_smoke(algorithm, config):
    sequence = "GGCACGAUCGGGCUCGCUGCCUUUUCGUCCGAGAGCUCGAA"
    results = algorithm.get_results(sequence, **config)

    for _, r in results.iterrows():
        assert isinstance(r.dot_bracket, str) and len(r.dot_bracket) == len(sequence)
        assert isinstance(r.energy, float)
        assert isinstance(r.stems, int)
