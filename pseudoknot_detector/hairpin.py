import argparse
import itertools

import ctypes
import pandas as pd

from pseudoknot_detector.grammars.hairpin import generate_grammar

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


def main():
    parser = argparse.ArgumentParser("rna_hairpin")
    parser.add_argument("sequence")
    parser.add_argument("--grammar", required=True)
    parser.add_argument("--allow-ug", default=False, action="store_true")
    parser.add_argument("--min-stems", default=MIN_HAIRPIN_STEMS, type=int)
    parser.add_argument("--min-size", default=MIN_HAIRPIN_SIZE, type=int)
    parser.add_argument("--max-per-loop", default=MAX_HAIRPINS_PER_LOOP, type=int)
    parser.add_argument("--max-bulge", default=MAX_HAIRPIN_BULGE, type=int)

    args = parser.parse_args()
    hairpin = HairpinDetector(
        grammar=args.grammar,
        allow_ug=args.allow_ug,
        min_stems=args.min_stems,
        min_size=args.min_size,
        max_per_loop=args.max_per_loop,
        max_bulge=args.max_bulge,
    )

    print("\n".join(hairpin.detect_hairpins(args.sequence)))


if __name__ == "__main__":
    main()
