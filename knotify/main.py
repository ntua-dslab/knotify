import argparse
from datetime import datetime

import pandas as pd

from knotify import rna_analysis
from knotify import hairpin
from knotify.criteria import apply_free_energy_and_stems_criterion
from knotify.energy.base import BaseEnergy
from knotify.energy.vienna import ViennaEnergy
from knotify.energy.pkenergy import PKEnergy
from knotify.parsers.base import BaseParser
from knotify.parsers.yaep import YaepParser
from knotify.parsers.bruteforce import BruteForceParser

from knotify.foreign import ipknot


def get_results(
    sequence: str,
    parser: BaseParser,
    allow_skip_final_au: bool = False,
    csv: str = None,
    max_stem_allow_smaller: int = 1,
    prune_early: bool = False,
    hairpin_grammar: str = None,
    hairpin_allow_ug: bool = False,
    min_hairpin_size: int = hairpin.MIN_HAIRPIN_SIZE,
    min_hairpin_stems: int = hairpin.MIN_HAIRPIN_STEMS,
    max_hairpins_per_loop: int = hairpin.MAX_HAIRPINS_PER_LOOP,
    max_hairpin_bulge: int = hairpin.MAX_HAIRPIN_BULGE,
    energy: BaseEnergy = ViennaEnergy(),
    *args,
    **kwargs,
) -> pd.DataFrame:
    """
    Analyze RNA sequence, and predict structure. Return data frame of results
    """
    sequence = sequence.lower()
    knot_dict_list = rna_analysis.StringAnalyser(
        input_string=sequence,
        parser=parser,
    ).get_pseudoknots(
        max_stem_allow_smaller=max_stem_allow_smaller,
        prune_early=prune_early,
        allow_skip_final_au=allow_skip_final_au,
    )

    data = pd.DataFrame(knot_dict_list)
    if csv is not None:
        data.to_csv(csv)

    data = apply_free_energy_and_stems_criterion(
        data,
        sequence,
        max_stem_allow_smaller=max_stem_allow_smaller,
        energy=energy,
    )

    if hairpin_grammar is None:
        return data

    data = hairpin.find_hairpins(
        sequence,
        data,
        hairpin_grammar,
        allow_ug=hairpin_allow_ug,
        min_size=min_hairpin_size,
        min_stems=min_hairpin_stems,
        max_bulge=max_hairpin_bulge,
        max_per_loop=max_hairpins_per_loop,
    )
    data = apply_free_energy_and_stems_criterion(
        data,
        sequence,
        max_stem_allow_smaller=max_stem_allow_smaller,
        energy=energy,
    )

    return data


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
    parser.add_argument("--max-loop-size", default=100, type=int)
    parser.add_argument("--min-loop-size", default=1, type=int)
    parser.add_argument("--max-stem-allow-smaller", default=2, type=int)
    parser.add_argument("--prune-early", default=False, action="store_true")

    # hairpin arguments
    parser.add_argument("--hairpin-grammar")
    parser.add_argument("--hairpin-allow-ug", default=None, action="store_true")
    parser.add_argument(
        "--min-hairpin-size", type=int, default=hairpin.MIN_HAIRPIN_SIZE
    )
    parser.add_argument(
        "--min-hairpin-stems", type=int, default=hairpin.MIN_HAIRPIN_STEMS
    )
    parser.add_argument(
        "--max-hairpins-per-loop", type=int, default=hairpin.MAX_HAIRPINS_PER_LOOP
    )
    parser.add_argument(
        "--max-hairpin-bulge", type=int, default=hairpin.MAX_HAIRPIN_BULGE
    )

    # energy arguments
    parser.add_argument("--energy", choices=["vienna", "pkenergy"], default="vienna")
    parser.add_argument("--pkenergy", default="./libpkenergy.so")
    parser.add_argument("--pkenergy-config-dir", default="pkenergy/hotknots/params")

    # overrides for other algorithms
    parser.add_argument("--algorithm", choices=["knotify", "ipknot"], default="knotify")
    parser.add_argument("--ipknot-executable", default="ipknot")

    return parser


def config_from_arguments(args: argparse.Namespace) -> dict:
    if args.hairpin_allow_ug is None:
        args.hairpin_allow_ug = args.allow_ug

    rna_parser_args = {
        "max_dd_size": args.max_dd_size,
        "min_dd_size": args.min_dd_size,
        "max_window_size": 2 * args.max_loop_size + 4,
        "min_window_size": 2 * args.min_loop_size + 4,
        "allow_ug": args.allow_ug,
    }
    if args.parser == "yaep":
        parser = YaepParser(args.yaep_library_path, **rna_parser_args)
    elif args.parser == "bruteforce":
        parser = BruteForceParser(args.bruteforce_library_path, **rna_parser_args)

    if args.energy == "vienna":
        energy = ViennaEnergy()
    elif args.energy == "pkenergy":
        energy = PKEnergy(args.pkenergy, args.pkenergy_config_dir)

    if args.algorithm == "knotify":
        algorithm = get_results
    elif args.algorithm == "ipknot":
        algorithm = ipknot.get_results

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
    }


def main():

    # define the argument parser
    parser = argument_parser()
    parser.add_argument("sequence")

    args = parser.parse_args()

    config = config_from_arguments(args)
    algorithm = config["algorithm"]

    start = datetime.now()
    results = algorithm(sequence=args.sequence.lower(), **config)
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
