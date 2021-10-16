import argparse

from knotify.parsers.ctypes import CTypesParser
from knotify.grammars.pseudoknot import generate_grammar


class YaepParser(CTypesParser):
    """
    Parser that uses a YAEP grammar to detect pseudoknots in a string.

    Reference C library implementation is in parsers/pseudoknot.c
    """

    def __init__(self, *args, **kwargs):
        super(YaepParser, self).__init__(*args, **kwargs)

    def get_options(self) -> str:
        return generate_grammar(self.allow_ug, self.max_dd_size)
