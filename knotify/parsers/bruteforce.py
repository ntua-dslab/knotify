from knotify.parsers.ctypes import CTypesParser
from knotify.grammars.pseudoknot import generate_grammar


class BruteForceParser(CTypesParser):
    """
    Parser that uses the pseudoknot grammar to detect pseudoknots in a string.
    """

    def __init__(self, *args, **kwargs):
        super(BruteForceParser, self).__init__(*args, **kwargs)

    def get_options(self) -> str:
        return ""
