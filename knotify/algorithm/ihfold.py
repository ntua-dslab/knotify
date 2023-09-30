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
import logging
import re
import subprocess

import pandas as pd

from knotify.algorithm.base import BaseAlgorithm


LOG = logging.getLogger(__name__)

PATTERN_V1 = re.compile(r"RES:\s([\.\{\}\(\)\[\]]*)\s+(.*)\n")
PATTERN_V2 = re.compile(r"Result_0:\s*([\.\{\}\(\)\[\]]*)\s+\((.*)\)\n")


def parse_output(stdout: str, pattern: re.Pattern) -> dict:
    """
    Parse output of HFold_iterative execution, return an object of the following form:

    {
        "dot_bracket": "...((((.[[[..))))..]]]...",
        "energy": 32.412,
        "stems": 0,
    }
    """
    dot_bracket, energy = pattern.findall(stdout)[0]
    return {
        "dot_bracket": dot_bracket,
        "stems": len(dot_bracket) - sum(x == "." for x in dot_bracket),
        "energy": float(energy),
    }


class IHFold(BaseAlgorithm):
    def get_results(self, sequence: str, ihfold_executable: str, *args, **kwargs):
        """
        Run HFold-Iterative algorithm

        {
            "dot_bracket": "...((((.[[[..))))..]]]...",
            "energy": 32.412,
            "stems": 0,
        }

        Example invocation of HFold_iterative:
        ```bash
        $ ./HFold_iterative \
            --s GGACGAGGAGCGCUGCAAGCGAGAGCCCAGGCUCGUCCGUUCAAACGGCGCUCA \
            --r ______________________________________________________
        method1
        method2
        method3
        method4
        Seq: GGACGAGGAGCGCUGCAAGCGAGAGCCCAGGCUCGUCCGUUCAAACGGCGCUCA
        RES: .......[[[[[[[[..[[[[.[[[[....]]]]...]]]]....]]]]]]]].  -14.94
        ```
        """
        cmd = [ihfold_executable, "--s", sequence.upper(), "--r", "_" * len(sequence)]
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if p.returncode != 0:
            LOG.warning(
                "ihfold invocation %s failed with return code %d", cmd, p.returncode
            )

        try:
            result = parse_output(p.stdout.decode(), PATTERN_V1)
        except IndexError:
            LOG.warning("ihfold invocation %s result could not be parsed", cmd)
            result = {"dot_bracket": "." * len(sequence), "stems": 0, "energy": 0}

        return pd.DataFrame([result])


class IHFoldV2(BaseAlgorithm):
    def get_results(self, sequence: str, ihfoldv2_executable: str, *args, **kwargs):
        """
        Run Iterative-HFold (ihfold-v2) algorithm

        Example invocation of Iterative-HFold:

        ```bash
        # ./Iterative-HFold CACCGUACCUAUUUAGGUUU
        ______((((____))))__
        Seq:          CACCGUACCUAUUUAGGUUU
        Restricted_0: ______((((____))))__
        Result_0:     ....[[((((]]..)))).. (-1.61)
        ```
        """

        cmd = [ihfoldv2_executable, sequence.upper()]
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if p.returncode != 0:
            LOG.warning(
                "ihfold invocation %s failed with return code %d", cmd, p.returncode
            )

        try:
            result = parse_output(p.stdout.decode(), PATTERN_V2)
        except IndexError:
            LOG.warning("ihfold invocation %s result could not be parsed", cmd)
            result = {"dot_bracket": "." * len(sequence), "stems": 0, "energy": 0}

        return pd.DataFrame([result])
