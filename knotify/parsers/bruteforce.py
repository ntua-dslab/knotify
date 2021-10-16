from knotify.parsers.ctypes import CTypesParser
from knotify.grammars.pseudoknot import generate_grammar


class BruteForceParser(CTypesParser):
    """
    Parser that uses brute force to detect all possible pseudoknot core stems.

    Reference C library implementation is in parsers/bruteforce.c
    """

    def __init__(self, *args, **kwargs):
        super(BruteForceParser, self).__init__(*args, **kwargs)

    def get_options(self) -> str:
        return ""
