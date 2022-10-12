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
import os

import pytest

from knotify.energy.pkenergy import PKEnergy

PKENERGY_SO = os.getenv("PKENERGY_SO", "./libpkenergy.so")
PKENERGY_PARAMS = os.getenv("PKENERGY_PARAMS", "pkenergy/hotknots/params")
PKENERGY_MODEL = os.getenv("PKENERGY_MODEL", "dp")


@pytest.mark.parametrize(
    "sequence, dot_bracket, model, result",
    [
        (
            "AAUGCAACUUUUAAAUAGUUUAUCUGUUAAGAUAAACCACCUAGGUUGCAUAUAUAAAAAAUAAAAGGUGCC",
            ".(((((((((...........................[[[[[)))))))))..............]]]]]..",
            "dp",
            -6.985999584197998,
        ),
        (
            "AAUGCAACUUUUAAAUAGUUUAUCUGUUAAGAUAAACCACCUAGGUUGCAUAUAUAAAAAAUAAAAGGUGCC",
            ".(((((((((...........................[[[[[)))))))))..............]]]]]..",
            "cc2006b",
            -11.74093246459961,
        ),
    ],
)
def test_pkenergy(sequence: str, dot_bracket: str, model: str, result: float):
    e = PKEnergy(PKENERGY_SO, PKENERGY_PARAMS, model)
    assert e.eval(sequence, dot_bracket) == result
