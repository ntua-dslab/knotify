#
# Copyright © 2021 Christos Pavlatos, George Rassias, Christos Andrikos,
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
import argparse
from datetime import datetime

from knotify import hairpin
from knotify.algorithm.ipknot import IPKnot
from knotify.algorithm.knotify import Knotify
from knotify.algorithm.knotty import Knotty
from knotify.algorithm.ihfold import IHFold
from knotify.algorithm.hotknots import HotKnots
from knotify.energy.vienna import ViennaEnergy
from knotify.energy.pkenergy import PKEnergy
from knotify.parsers.yaep import YaepParser
from knotify.parsers.bruteforce import BruteForceParser


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    # store intermediate results
    parser.add_argument("--csv")
    parser.add_argument("--results-csv")

    # pseudoknot arguments
    parser.add_argument("--parser", choices=["bruteforce", "yaep"])
    parser.add_argument("--yaep-library-path", default="./libpseudoknot.so")
    parser.add_argument("--bruteforce-library-path", default="./libbruteforce.so")
    parser.add_argument("--allow-ug", default=False, action="store_true")
    parser.add_argument("--allow-skip-final-au", default=False, action="store_true")
    parser.add_argument("--max-dd-size", default=2, type=int)
    parser.add_argument("--min-dd-size", default=0, type=int)
    parser.add_argument("--max-window-size", default=204, type=int)
    parser.add_argument("--min-window-size", default=6, type=int)
    parser.add_argument("--max-window-size-ratio", default=0, type=float)
    parser.add_argument("--min-window-size-ratio", default=0, type=float)
    parser.add_argument("--max-stem-allow-smaller", default=2, type=int)
    parser.add_argument("--prune-early", default=False, action="store_true")

    # hairpin arguments
    parser.add_argument("--hairpin-grammar")
    parser.add_argument("--hairpin-allow-ug", default=None, action="store_true")
    parser.add_argument("--min-hairpin-size", type=int, default=hairpin.MIN_HAIRPIN_SIZE)
    parser.add_argument("--min-hairpin-stems", type=int, default=hairpin.MIN_HAIRPIN_STEMS)
    parser.add_argument("--max-hairpins-per-loop", type=int, default=hairpin.MAX_HAIRPINS_PER_LOOP)
    parser.add_argument("--max-hairpin-bulge", type=int, default=hairpin.MAX_HAIRPIN_BULGE)

    # energy arguments
    parser.add_argument("--energy", choices=["vienna", "pkenergy"], default="vienna")
    parser.add_argument("--pkenergy", default="./libpkenergy.so")
    parser.add_argument("--pkenergy-config-dir", default="pkenergy/hotknots/params")

    # overrides for other algorithms
    parser.add_argument(
        "--algorithm",
        choices=["knotify", "ipknot", "knotty", "hotknots", "ihfold"],
        default="knotify",
    )
    parser.add_argument("--ipknot-executable", default="./.ipknot/ipknot")
    parser.add_argument("--knotty-executable", default="./.knotty/knotty")
    parser.add_argument("--ihfold-executable", default="./.iterative-hfold/HFold_iterative")
    parser.add_argument("--hotknots-dir", default="./.hotknots/HotKnots_v2.0")

    return parser


def config_from_arguments(args: argparse.Namespace) -> dict:
    if args.hairpin_allow_ug is None:
        args.hairpin_allow_ug = args.allow_ug

    rna_parser_args = {
        "max_dd_size": args.max_dd_size,
        "min_dd_size": args.min_dd_size,
        "max_window_size": args.max_window_size,
        "min_window_size": args.min_window_size,
        "max_window_size_ratio": args.min_window_size_ratio,
        "min_window_size_ratio": args.max_window_size_ratio,
        "allow_ug": args.allow_ug,
    }
    parser = None
    if args.parser == "yaep":
        parser = YaepParser(args.yaep_library_path, **rna_parser_args)
    elif args.parser == "bruteforce":
        parser = BruteForceParser(args.bruteforce_library_path, **rna_parser_args)

    energy = None
    if args.energy == "vienna":
        energy = ViennaEnergy()
    elif args.energy == "pkenergy":
        energy = PKEnergy(args.pkenergy, args.pkenergy_config_dir)

    algorithm = None
    if args.algorithm == "knotify":
        algorithm = Knotify()
    elif args.algorithm == "ipknot":
        algorithm = IPKnot()
    elif args.algorithm == "knotty":
        algorithm = Knotty()
    elif args.algorithm == "ihfold":
        algorithm = IHFold()
    elif args.algorithm == "hotknots":
        algorithm = HotKnots()

    return {
        "algorithm": algorithm,
        "parser": parser,
        "csv": args.csv,
        "allow_ug": args.allow_ug,
        "allow_skip_final_au": args.allow_skip_final_au,
        "max_stem_allow_smaller": args.max_stem_allow_smaller,
        "prune_early": args.prune_early,
        "hairpin_grammar": args.hairpin_grammar,
        "hairpin_allow_ug": args.hairpin_allow_ug,
        "min_hairpin_size": args.min_hairpin_size,
        "min_hairpin_stems": args.min_hairpin_stems,
        "max_hairpin_bulge": args.max_hairpin_bulge,
        "max_hairpins_per_loop": args.max_hairpins_per_loop,
        "energy": energy,
        "ipknot_executable": args.ipknot_executable,
        "knotty_executable": args.knotty_executable,
        "ihfold_executable": args.ihfold_executable,
        "hotknots_dir": args.hotknots_dir,
    }


def main():

    # define the argument parser
    parser = argument_parser()
    parser.add_argument("sequence")

    args = parser.parse_args()

    config = config_from_arguments(args)
    algorithm = config["algorithm"]

    start = datetime.now()
    results = algorithm.get_results(sequence=args.sequence.lower(), **config)
    duration = datetime.now() - start

    if args.results_csv:
        results[["dot_bracket", "stems", "energy"]].to_csv(args.results_csv, index=None)

    chosen = results.loc[0]

    print("Sequence: ", args.sequence)
    print("Structure:", chosen.dot_bracket)
    print("Energy:", chosen.energy)
    print("Duration:", duration.total_seconds(), "s")


if __name__ == "__main__":
    main()
