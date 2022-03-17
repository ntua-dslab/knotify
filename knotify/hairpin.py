#
# Copyright © 2022 Christos Pavlatos, George Rassias, Christos Andrikos,
#                  Evangelos Makris, Aggelos Kolaitis
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the “Software”), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
import itertools

import ctypes
import pandas as pd
from oslo_config import cfg

from knotify.grammars.hairpin import generate_grammar

MIN_HAIRPIN_STEMS = 3
MIN_HAIRPIN_SIZE = 3
MAX_HAIRPIN_BULGE = 0
MAX_HAIRPINS_PER_LOOP = 1


class HairpinDetector:
    """
    Load detector from a dynamic library. The detector function should be defined
    in C code as:

    char *detect_hairpins(
        char* grammar,
        char *sequence,
        int min_stems,
        int min_size,
        int max_per_loop,
        int max_bulge
    )
    """

    def __init__(
        self,
        grammar: str,
        allow_ug: bool,
        min_stems: int = MIN_HAIRPIN_STEMS,
        min_size: int = MIN_HAIRPIN_SIZE,
        max_per_loop: bool = False,
        max_bulge: int = MAX_HAIRPIN_BULGE,
    ):
        self.grammar = grammar
        self.definition = generate_grammar(allow_ug=allow_ug)
        self.min_stems = min_stems
        self.min_size = min_size
        self.max_per_loop = max_per_loop
        self.max_bulge = max_bulge

    def detect_hairpins(self, loop_sequence: str):
        if not loop_sequence:
            return [""]

        if not hasattr(self, "lib"):
            self.lib = ctypes.CDLL(self.grammar)
            self.lib.detect_hairpins.restype = ctypes.c_char_p

        loops = (
            self.lib.detect_hairpins(
                ctypes.c_char_p(self.definition.encode()),
                ctypes.c_char_p(loop_sequence.lower().encode()),
                self.min_stems,
                self.min_size,
                self.max_per_loop,
                self.max_bulge,
            )
            .decode()
            .rstrip("\n")
            .split("\n")
        )

        # drop duplicates
        return list(set(loops))


def get_loop_indices(dot_bracket: str):
    """
    Return indices for left and right loop sequences.
    """
    return (
        dot_bracket.rindex("(") + 1,
        dot_bracket.index("["),
        dot_bracket.rindex(")") + 1,
        dot_bracket.index("]"),
    )


def find_hairpins_in_loop(detector: HairpinDetector, loop_sequence: str) -> set:
    """
    Accepts an RNA loop sequence.

    Returns a set of loops including hairpins. This function is responsible
    for implementing and enforcing any pruning criteria on the possible
    hairpin loops.

    It makes sure that a set containing only the empty string is returned
    if the loop size is zero. It also makes sure that a loop with no
    hairpins is always included in the result.
    """
    if loop_sequence == "":
        return set([""])

    all_loops = detector.detect_hairpins(loop_sequence)

    # TODO(akolaitis): decide and apply other limiting criteria
    loops = set(all_loops + ["." * len(loop_sequence)])
    if len(loops) > 1 and "" in loops:
        loops.remove("")

    return loops


def dot_bracket_to_record(dot_bracket: str) -> tuple:
    """
    Calculate left loop stems, right loop stems and dd of dot bracket.

    This assumes that the dot bracket is valid contains core stems only.
    """
    left_loop_stems = sum(x == "(" for x in dot_bracket) - 1
    right_loop_stems = sum(x == "[" for x in dot_bracket) - 1
    dd = dot_bracket.index(")") - dot_bracket.rindex("[") - 1

    return left_loop_stems, right_loop_stems, dd


def find_hairpins_in_dot_bracket(
    detector: HairpinDetector, sequence: str, dot_bracket: str
) -> pd.DataFrame:
    """
    Accept an RNA sequence and a possible dot bracket. This function
    finds the left and right loop indices, gets a list of all possible
    hairpin combinations and returns a pandas.DataFrame with records
    of the following format:

    {
        "dot_bracket": "....(((([[[))))...]]]",
        "left_loop_stems": <core left loop stems>,
        "right_loop_stems": <core right loop stems>,
        "dd": <dd length>,
    }
    """
    lstart, lend, rstart, rend = get_loop_indices(dot_bracket)
    left_loops = find_hairpins_in_loop(detector, sequence[lstart:lend])
    right_loops = find_hairpins_in_loop(detector, sequence[rstart:rend])

    left_loop_stems, right_loop_stems, dd = dot_bracket_to_record(dot_bracket)

    return pd.DataFrame.from_records(
        [
            {
                "dot_bracket": (
                    dot_bracket[:lstart]
                    + left_loop
                    + dot_bracket[lend:rstart]
                    + right_loop
                    + dot_bracket[rend:]
                ),
                "left_loop_stems": left_loop_stems,
                "right_loop_stems": right_loop_stems,
                "dd": dd,
            }
            for left_loop, right_loop in itertools.product(left_loops, right_loops)
        ]
    )


def find_hairpins(
    sequence: str,
    data: pd.DataFrame,
    hairpin_grammar: str,
    allow_ug: bool = False,
    min_stems: int = MIN_HAIRPIN_STEMS,
    min_size: int = MIN_HAIRPIN_SIZE,
    max_per_loop: int = MAX_HAIRPINS_PER_LOOP,
    max_bulge: int = MAX_HAIRPIN_BULGE,
) -> pd.DataFrame:
    """
    For each row in the specified data frame, try to find hairpins in each loop.
    Generate a new data frame with all possible combinations.
    """
    detector = HairpinDetector(
        grammar=hairpin_grammar,
        allow_ug=allow_ug,
        min_stems=min_stems,
        min_size=min_size,
        max_per_loop=max_per_loop,
        max_bulge=max_bulge,
    )

    parts = [
        find_hairpins_in_dot_bracket(detector, sequence, row.dot_bracket)
        for _, row in data.iterrows()
    ]

    return pd.concat(parts, ignore_index=True)


OPTS = [
    cfg.StrOpt("sequence"),
    cfg.StrOpt("grammar", default="./libhairpin.so"),
    cfg.BoolOpt("allow-ug", default=False),
    cfg.IntOpt("min-stems", default=MIN_HAIRPIN_STEMS),
    cfg.IntOpt("min-size", default=MIN_HAIRPIN_SIZE),
    cfg.IntOpt("max-per-loop", default=MAX_HAIRPINS_PER_LOOP),
    cfg.IntOpt("max-bulge", default=MAX_HAIRPIN_BULGE),
]


def main():
    options = cfg.ConfigOpts()
    options.register_cli_opts(OPTS)
    options()

    if not options.sequence:
        print("Missing required argument --sequence")
        sys.exit(1)

    hairpin = HairpinDetector(
        grammar=options.grammar,
        allow_ug=options.allow_ug,
        min_stems=options.min_stems,
        min_size=options.min_size,
        max_per_loop=options.max_per_loop,
        max_bulge=options.max_bulge,
    )

    print("\n".join(hairpin.detect_hairpins(options.sequence)))


if __name__ == "__main__":
    main()
