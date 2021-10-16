import os

import pytest

from knotify.parsers.bruteforce import BruteForceParser
from knotify.parsers.yaep import YaepParser

PSEUDOKNOT = os.getenv("PSEUDOKNOT_SO", "./libpseudoknot.so")
BRUTEFORCE = os.getenv("BRUTEFORCE_SO", "./libbruteforce.so")


def for_each_parser(variable_names):
    def wrapped(f):
        return pytest.mark.parametrize(
            variable_names,
            [
                (YaepParser, PSEUDOKNOT),
                (BruteForceParser, BRUTEFORCE),
            ],
        )(f)

    return wrapped
