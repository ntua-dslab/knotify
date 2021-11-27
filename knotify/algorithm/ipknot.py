import logging
import re
import subprocess
import tempfile

import pandas as pd

from knotify.algorithm.base import BaseAlgorithm


LOG = logging.getLogger(__name__)

ENERGY_REGEX = re.compile(r"^> case\s*\(e=(.*)\)$")


def normalize_dot_bracket(bracket: str) -> str:
    """
    Return a normalized version of the dot bracket. First pseudoknot pair
    always starts with a parenthesis.

    ```
    normalize_dot_bracket("..((..[[[..)).]]].") == "..((..[[[..)).]]]."
    normalize_dot_bracket("..[[..(((..]].))).") == "..((..[[[..)).]]]."
    ```
    """
    if "[" not in bracket or "(" not in bracket:
        return bracket
    if bracket.index("[") > bracket.index("("):
        return bracket

    return (
        bracket.replace("[", "A")
        .replace("]", "B")
        .replace("(", "[")
        .replace(")", "]")
        .replace("A", "(")
        .replace("B", ")")
    )


def parse_ipknot_output(stdout: str) -> dict:
    """
    Parse output of ipknot execution, return an object of the following form:

    {
        "dot_bracket": "...((((.[[[..))))..]]]...",
        "energy": 32.412,
        "stems": 0,
    }
    """
    lines = stdout.strip().split("\n")
    try:
        energy = float(ENERGY_REGEX.findall(lines[0])[0])
    except (IndexError, ValueError, TypeError) as e:
        energy = 1000
        LOG.warning("failed to parse energy from %s", lines)

    dot_bracket = normalize_dot_bracket(lines[-1])
    return {
        "dot_bracket": dot_bracket,
        "stems": len(dot_bracket) - sum(x == "." for x in dot_bracket),
        "energy": energy,
    }


class IPKnot(BaseAlgorithm):
    def get_results(self, sequence: str, ipknot_executable: str, *args, **kwargs):
        """
        Run ipknot algorithm

        {
            "dot_bracket": "...((((.[[[..))))..]]]...",
            "energy": 32.412,
            "stems": 0,
        }

        Example invocation of ipknot:
        ```bash
        $ cat file
        > case
        AAAAAACUAAUAGAGGGGGGACUUAGCGCCCCCCAAACCGUAACCCC

        $ ./ipknot -E file
        > case (e=-26.934)
        AAAAAACUAAUAGAGGGGGGACUUAGCGCCCCCCAAACCGUAACCCC
        ..............((((((........)))))).............
        ```
        """
        with tempfile.NamedTemporaryFile() as f:
            f.write("> case\n{}".format(sequence.upper()).encode())
            f.flush()
            cmd = [ipknot_executable, "-E", f.name]
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if p.returncode != 0:
            LOG.warning(
                "ipknot invocation %s for %s failed with return code %d",
                cmd,
                sequence,
                p.returncode,
            )

        return pd.DataFrame([parse_ipknot_output(p.stdout.decode())])
