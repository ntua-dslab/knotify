import os
from typing import List

import pytest

from pseudoknot_detector.rna_parser import PseudoknotDetector

PSEUDOKNOT = os.getenv("PSEUDOKNOT_SO", "./libpseudoknot.so")


@pytest.mark.parametrize(
    "dd_size, sequence, result",
    [
        (1, "caaaaaggaaaaac", ["5, 0"]),
        (1, "caaaaagagaaaaac", ["5, 1"]),
        (1, "caaaaagaagaaaaac", [""]),
        (2, "caaaaagaagaaaaac", ["5, 2"]),
        (2, "caaaaagaaagaaaaac", [""]),
        (3, "caaaaagaaagaaaaac", ["5, 3"]),
    ],
)
def test_max_dd_size(dd_size: int, sequence: str, result: List[str]):
    parser = PseudoknotDetector(grammar=PSEUDOKNOT, max_dd_size=dd_size, allow_ug=False)
    assert parser.detect_pseudoknots(sequence) == result
