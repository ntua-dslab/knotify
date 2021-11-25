import logging
import re
import subprocess
import tempfile

import pandas as pd


LOG = logging.getLogger(__name__)

ENERGY_REGEX = re.compile(r"^> case\s*\(e=(.*)\)$")


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

    return {
        "dot_bracket": lines[-1],
        "stems": -1,  # TODO: count '([])' occurences
        "energy": energy,
    }


def get_results(sequence: str, ipknot_executable: str, *args, **kwargs):
    """
    Parse output of ipknot for an RNA sequence, return the result in a Pandas
    dataframe with rows of the following format:

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
    with tempfile.TemporaryNamedFile() as f:
        f.write(b"> case\n{}".format(sequence))

    cmd = [ipknot_executable, "-E", f.name]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if p.returncode != 0:
        LOG.warning(
            "ipknot invocation %s failed with return code %d", cmd, p.returncode
        )

    return pd.DataFrame([parse_ipknot_output(p.stdout.decode())])
