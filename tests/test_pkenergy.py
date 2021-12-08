#
# Copyright © 2021 Christos Pavlatos, George Rassias, Christos Andrikos,
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
    assert e.eval(sequence, dot_bracket) == result
