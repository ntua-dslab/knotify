import argparse

from knotify.parsers.ctypes import CTypesParser
from knotify.grammars.pseudoknot import generate_grammar


class YaepParser(CTypesParser):
    """
    Parser that uses the pseudoknot grammar to detect pseudoknots in a string.
    """

    def __init__(self, *args, **kwargs):
        super(YaepParser, self).__init__(*args, **kwargs)

    def get_options(self) -> str:
        return generate_grammar(self.allow_ug, self.max_dd_size)


def main():
    prog = argparse.ArgumentParser("rna_parser")
    prog.add_argument("sequence", type=str)
    prog.add_argument("--library-path", type=str, required=True)
    prog.add_argument("--max-dd-size", type=int, default=2)
    prog.add_argument("--min-dd-size", type=int, default=0)
    prog.add_argument("--allow-ug", action="store_true", default=False)

    args = prog.parse_args()
    parser = YaepParser(
        library_path=args.library_path,
        max_dd_size=args.max_dd_size,
        allow_ug=args.allow_ug,
        min_dd_size=args.min_dd_size,
    )

    print(parser.detect_pseudoknots(args.sequence))
