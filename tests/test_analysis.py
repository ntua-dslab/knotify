import pytest

from knotify.rna_analysis import StringAnalyser
from tests.utils import for_each_parser


@pytest.mark.parametrize(
    "sequence, brackets",
    [
        (
            "AUCCUUUUCA",
            [
                "(.....[).]",
                "(...[..).]",
                "(....[)..]",
                "(...[)...]",
                "(....[.).]",
                "(...[.)..]",
            ],
        ),
        (
            "AGAUUGACU",
            [
                "(..[).]..",
                "(.[.)...]",
                "(.[)....]",
            ],
        ),
    ],
)
@for_each_parser("parser, library_path")
def test_get_pseudoknots(parser, library_path, sequence, brackets):
    s = StringAnalyser(
        sequence, parser=parser(library_path=library_path)
    ).get_pseudoknots()

    assert len(s) == len(brackets)
    for a in s:
        assert a["dot_bracket"] in brackets
