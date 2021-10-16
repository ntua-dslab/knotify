import argparse
import subprocess
import os

import ctypes

from knotify.grammars.pseudoknot import generate_grammar


class PseudoknotDetector:
    """
    Load parser from a dynamic library. The parser function should be defined
    in C code as

    ```
    char *detect_pseudoknots(
        char* grammar, char *sequence, int max_dd_size, int min_dd_size
    );
    ```
    """

    def __init__(
        self,
        grammar: str,
        max_dd_size: int,
        allow_ug: bool,
        min_dd_size: int = 0,
        max_window_size: int = 100,
        min_window_size: int = 6,
    ):
        self.grammar = grammar
        self.max_dd_size = max_dd_size
        self.min_dd_size = min_dd_size
        self.max_window_size = max_window_size
        self.min_window_size = min_window_size
        self.definition = generate_grammar(allow_ug=allow_ug, max_dd_size=max_dd_size)

    def detect_pseudoknots(self, x):
        if not hasattr(self, "lib"):
            self.lib = ctypes.CDLL(self.grammar)
            self.lib.detect_pseudoknots.restype = ctypes.c_char_p

        # TODO(akolaitis): currently lists of results are returned as a C string,
        # and then split and parsed as integers with slow Python code.
        #
        # Consider changing the definition of the C function to return an array
        # of integers instead.
        #
        # This would require updating the implementation of StringAnalyzer,
        # specifically the `.dd()` and `.left_loop_boundary()` methods.
        return (
            self.lib.detect_pseudoknots(
                ctypes.c_char_p(self.definition.encode()),
                ctypes.c_char_p(x.lower().encode()),
                ctypes.c_int(self.max_dd_size),
                ctypes.c_int(self.min_dd_size),
                ctypes.c_int(self.max_window_size),
                ctypes.c_int(self.min_window_size),
            )
            .decode()
            .rstrip("\n")
            .split("\n")
        )


def main():
    prog = argparse.ArgumentParser("rna_parser")
    prog.add_argument("sequence", type=str)
    prog.add_argument("--grammar", type=str, required=True)
    prog.add_argument("--max-dd-size", type=int, default=2)
    prog.add_argument("--min-dd-size", type=int, default=0)
    prog.add_argument("--allow-ug", action="store_true", default=False)

    args = prog.parse_args()
    parser = PseudoknotDetector(
        args.grammar,
        max_dd_size=args.max_dd_size,
        allow_ug=args.allow_ug,
        min_dd_size=args.min_dd_size,
    )

    print(parser.detect_pseudoknots(args.sequence))
