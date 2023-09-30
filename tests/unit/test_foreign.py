#
# Copyright © 2022 Christos Pavlatos, George Rassias, Christos Andrikos,
#                  Evangelos Makris, Aggelos Kolaitis
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the “Software”), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
import pathlib
import pytest

from knotify.algorithm import ipknot
from knotify.algorithm import knotty
from knotify.algorithm import ihfold
from knotify.algorithm import hotknots

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


IHFOLD_OUTPUT = "Seq: AUGAAAG\nRES: .([.)].  -3.2\n"
IHFOLD_OUTPUT_ALL = "Seq: AUGAAAAG\nRES: .([{})].  -3.21\n"


@pytest.mark.parametrize(
    "stdout, result",
    [
        (IHFOLD_OUTPUT, {"dot_bracket": ".([.)].", "stems": 4, "energy": -3.2}),
        (IHFOLD_OUTPUT_ALL, {"dot_bracket": ".([{})].", "stems": 6, "energy": -3.21}),
    ],
)
def test_ihfold_parse(stdout, result):
    assert ihfold.parse_output(stdout, ihfold.PATTERN_V1) == result


IHFOLDV2_OUTPUT = """
______((((____))))__
Seq:          CACCGUACCUAUUUAGGUUU
Restricted_0: ______((((____))))__
Result_0:     ....[[((((]]..)))).. (-1.61)
"""


@pytest.mark.parametrize(
    "stdout, result",
    [
        (
            IHFOLDV2_OUTPUT,
            {"dot_bracket": "....[[((((]]..))))..", "stems": 12, "energy": -1.61},
        ),
    ],
)
def test_ihfoldv2_parse(stdout, result):
    assert ihfold.parse_output(stdout, ihfold.PATTERN_V2) == result


HOTKNOTS_OUTPUT = (
    "Total number of RNA structures: 1\nSeq: AUGAAAG\nS0:  .([.)].\t-3.2\n"
)
HOTKNOTS_OUTPUT_ALL = (
    "Total number of RNA structures: 1\nSeq: AUGAAAG\nS0:  .([{})].\t-3.21\n"
)
HOTKNOTS_HARD = """Total number of RNA structures: 20
Seq: GGCACGAUCGGGCUCGCUGCCUUUUCGUCCGAGAGCUCGAA
S0:  ((([[[[{{{{{{{{)))......]]]]....}}}}}}}}.	-10.16
S1:  ...(((([[[[[[[[.........))))....]]]]]]]].	-8.73
S2:  .......((((((((...(.(.....).)...)))))))).	-8.44
S3:  ((((........[[[[.))))........]]]]........	-7.97
S4:  .......((((((((.................)))))))).	-7.84
S5:  ((([[[[..[[[[..)))]]]]..]]]].(((....)))..	-7.61
S6:  ((((...[[[[[[[[..))))...........]]]]]]]].	-7.45
S7:  ((((....[[[[{{{{.))))..]]]]..}}}}........	-7.33
S8:  (((....((((((.............))))))..)))....	-6.89
S9:  .......((((((..[[[........)))))).]]].....	-6.55
S10:  (((((((......))).))))..((((..........))))	-6.47
S11:  (((((((..((((.....))))..))))......)))....	-6.31
S12:  ...((((..((((.....))))..)))).(((....)))..	-6.28
S13:  ((((........[[[[.))))..((((..]]]]....))))	-6.22
S14:  ((((....[[[[.....))))..]]]]..(((....)))..	-5.65
S15:  .......((((((.............)))))).........	-5.58
S16:  (((......[[[[..)))]]]].((((..........))))	-5.00
S17:  ...((((..((((..[[[))))..)))).....]]].....	-4.66
S18:  ((((...[[[[[[....)))).....]]]]]].........	-4.59
S19:  .........((((..[[[)))).((((......]]].))))	-4.55
"""


@pytest.mark.parametrize(
    "stdout, result",
    [
        (
            HOTKNOTS_HARD,
            {
                "dot_bracket": "((([[[[{{{{{{{{)))......]]]]....}}}}}}}}.",
                "stems": 30,
                "energy": -10.16,
            },
        ),
        (HOTKNOTS_OUTPUT, {"dot_bracket": ".([.)].", "stems": 4, "energy": -3.2}),
        (HOTKNOTS_OUTPUT_ALL, {"dot_bracket": ".([{})].", "stems": 6, "energy": -3.21}),
    ],
)
def test_hotknots_parse(stdout, result):
    assert hotknots.parse_hotknots_output(stdout) == result


@pytest.mark.parametrize(
    "algorithm, executable, path",
    [
        (knotty.Knotty(), "knotty_executable", "./.knotty/knotty"),
        (ipknot.IPKnot(), "ipknot_executable", "./.ipknot/ipknot"),
        (hotknots.HotKnots(), "hotknots_dir", "./.hotknots/HotKnots_v2.0"),
        (ihfold.IHFold(), "ihfold_executable", "./.ihfold/v1/HFold_iterative"),
        (ihfold.IHFoldV2(), "ihfoldv2_executable", "./.ihfold/v2/Iterative-HFold"),
    ],
)
def test_foreign_smoke(algorithm, executable, path):
    if not pathlib.Path(path).exists():
        pytest.skip(f"Missing binary {path}")

    sequence = "GGCACGAUCGGGCUCGCUGCCUUUUCGUCCGAGAGCUCGAA"
    results = algorithm.get_results(sequence, **{executable: path})

    for _, r in results.iterrows():
        assert isinstance(r.dot_bracket, str) and len(r.dot_bracket) == len(sequence)
        assert isinstance(r.energy, float)
        assert isinstance(r.stems, int)
