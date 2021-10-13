import argparse
from datetime import datetime, timedelta
from typing import Tuple

import pandas as pd

from knotify import rna_analysis
from knotify import hairpin
from knotify.energy import apply_free_energy_and_stems_criterion


def get_results(
    sequence: str,
    grammar: str = None,
    max_dd_size: int = 2,
    min_dd_size: int = 0,
    allow_ug: bool = False,
    allow_skip_final_au: bool = False,
    max_loop_size: int = rna_analysis.MAX_LOOP_SIZE,
    save_csv: str = None,
    max_stem_allow_smaller: int = 1,
    prune_early: bool = False,
    hairpin_grammar: str = None,
    min_hairpin_size: int = hairpin.MIN_HAIRPIN_SIZE,
    min_hairpin_stems: int = hairpin.MIN_HAIRPIN_STEMS,
    max_hairpins_per_loop: int = hairpin.MAX_HAIRPINS_PER_LOOP,
    max_hairpin_bulge: int = hairpin.MAX_HAIRPIN_BULGE,
    energy_eval: callable = ViennaEnergy().energy_eval,
) -> pd.DataFrame:
    """
    Analyze RNA sequence, and predict structure. Return data frame of results
    """
    sequence = sequence.lower()
    knot_dict_list = (
        rna_analysis.StringAnalyser(
            input_string=sequence,
            grammar=grammar,
            max_loop_size=max_loop_size,
            max_dd_size=max_dd_size,
            min_dd_size=min_dd_size,
            allow_ug=allow_ug,
        )
        .get_window_boundaries()
        .generate_trees_in_parallel()
        .get_pseudoknots(
            max_stem_allow_smaller=max_stem_allow_smaller,
            prune_early=prune_early,
            allow_skip_final_au=allow_skip_final_au,
        )
    )

    data = pd.DataFrame(knot_dict_list)
    if save_csv is not None:
        data.to_csv(save_csv)

    data = apply_free_energy_and_stems_criterion(
        data,
        sequence,
        max_stem_allow_smaller=max_stem_allow_smaller,
        energy_eval=energy_eval,
    )

    if hairpin_grammar is None:
        return data

    data = hairpin.find_hairpins(
        sequence,
        data,
        hairpin_grammar,
        allow_ug=allow_ug,
        min_size=min_hairpin_size,
        min_stems=min_hairpin_stems,
        max_bulge=max_hairpin_bulge,
        max_per_loop=max_hairpins_per_loop,
    )
    data = apply_free_energy_and_stems_criterion(
        data,
        sequence,
        max_stem_allow_smaller=max_stem_allow_smaller,
        energy_eval=energy_eval,
    )

    return data


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    # store intermediate results
    parser.add_argument("--csv")
    parser.add_argument("--results-csv")

    # pseudoknot arguments
    parser.add_argument("--grammar")
    parser.add_argument("--allow-ug", default=False, action="store_true")
    parser.add_argument("--allow-skip-final-au", default=False, action="store_true")
    parser.add_argument("--max-dd-size", default=2, type=int)
    parser.add_argument("--min-dd-size", default=0, type=int)
    parser.add_argument("--max-loop-size", default=100, type=int)
    parser.add_argument("--max-stem-allow-smaller", default=2, type=int)
    parser.add_argument("--prune-early", default=False, action="store_true")

    # hairpin arguments
    parser.add_argument("--hairpin-grammar")
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
    return parser


def main():

    # define the argument parser
    parser = argument_parser()
    parser.add_argument("sequence")

    args = parser.parse_args()

    start = datetime.now()
    results = get_results(
        sequence=args.sequence.lower(),
        grammar=args.grammar,
        save_csv=args.csv,
        allow_ug=args.allow_ug,
        allow_skip_final_au=args.allow_skip_final_au,
        max_dd_size=args.max_dd_size,
        min_dd_size=args.min_dd_size,
        max_loop_size=args.max_loop_size,
        max_stem_allow_smaller=args.max_stem_allow_smaller,
        prune_early=args.prune_early,
        hairpin_grammar=args.hairpin_grammar,
        min_hairpin_size=args.min_hairpin_size,
        min_hairpin_stems=args.min_hairpin_stems,
        max_hairpin_bulge=args.max_hairpin_bulge,
        max_hairpins_per_loop=args.max_hairpins_per_loop,
    )
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
