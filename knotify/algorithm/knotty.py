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

RESULT_REGEX = re.compile(r"RES:\s([\.\{\}\(\)\[\]]*)\s+(.*)$")


def parse_knotty_output(stdout: str) -> dict:
    """
    Parse output of ipknot execution, return an object of the following form:

    {
        "dot_bracket": "...((((.[[[..))))..]]]...",
        "energy": 32.412,
        "stems": 0,
    }
    """
    dot_bracket, energy = RESULT_REGEX.findall(stdout)[0]
    return {
        "dot_bracket": dot_bracket,
        "stems": len(dot_bracket) - sum(x == "." for x in dot_bracket),
        "energy": float(energy),
    }


class Knotty(BaseAlgorithm):
    def get_results(self, sequence: str, knotty_executable: str, *args, **kwargs):
        """
        Run knotty algorithm

        {
            "dot_bracket": "...((((.[[[..))))..]]]...",
            "energy": 32.412,
            "stems": 0,
        }

        Example invocation of knotty:
        ```bash
        $ ./Knotty/knotty GGACGAGGAGCGCUGCAAGCGAGAGCCCAGGCUCGUCCGUUCAAACGGCGCUCA
        Seq: GGACGAGGAGCGCUGCAAGCGAGAGCCCAGGCUCGUCCGUUCAAACGGCGCUCA
        RES: ((((((([[[[[[[[.......[[[[.....)))))))]]]]...]]]]]]]].  -14.78
        ```
        """
        cmd = [knotty_executable, sequence]
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if p.returncode != 0:
            LOG.warning(
                "knotty invocation %s failed with return code %d", cmd, p.returncode
            )

        return pd.DataFrame([parse_knotty_output(p.stdout.decode())])
