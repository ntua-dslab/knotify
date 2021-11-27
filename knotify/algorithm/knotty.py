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
