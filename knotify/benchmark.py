import json
import sys
import logging

import yaml

from knotify.main import get_results, argument_parser
from knotify import scoring


LOG = logging.getLogger(__name__)


def main():
    parser = argument_parser()
    parser.add_argument("--cases", required=True)
    parser.add_argument("--only", type=int, nargs="*")
    parser.add_argument("--correct-stems-slack", type=int, default=0)
    parser.add_argument("--verbose", action="store_true", default=False)
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

        LOG.debug("Command line: %s", sys.argv)
        LOG.debug("Run configuration: %s", args.__dict__)

    with open(args.cases, "r") as fin:
        cases = yaml.safe_load(fin.read())

    only = args.only or list(range(0, len(cases)))

    out = {
        "results": [],
        "totals": {
            "new_correct": 0,
            "correct": 0,
            "same_prediction": 0,
            "correct_core_stems": 0,
            "truth_in_candidates": 0,
            "count": len(only),
        },
    }

    for loop, idx in enumerate(only):
        case = cases[idx]
        results = get_results(
            case["case"].lower(),
            grammar=args.grammar,
            save_csv=args.csv,
            allow_ug=args.allow_ug,
            max_dd_size=args.max_dd_size,
            max_loop_size=args.max_loop_size,
            max_stem_allow_smaller=args.max_stem_allow_smaller,
            prune_early=args.prune_early,
        )
        predictions = results[["dot_bracket", "energy", "stems"]].to_dict(
            orient="records"
        )
        dot_bracket = predictions[0]["dot_bracket"]

        new_correct = case["truth"] == dot_bracket
        if "pred" in case:
            correct = case["pred"] == case["truth"]
            same_prediction = case["pred"] == dot_bracket
        else:
            correct, same_prediction = 0, 0

        correct_core_stems = scoring.get_correct_core_stems(
            case["truth"], dot_bracket, slack=args.correct_stems_slack
        )
        confusion_matrix = scoring.get_confusion_matrix(
            case["truth"], dot_bracket, slack=0
        )
        truth_in_candidates = case["truth"] in [p["dot_bracket"] for p in predictions]

        out["totals"]["new_correct"] += new_correct
        out["totals"]["correct"] += correct
        out["totals"]["same_prediction"] += same_prediction
        out["totals"]["correct_core_stems"] += correct_core_stems == 2
        out["totals"]["truth_in_candidates"] += truth_in_candidates

        LOG.info(
            "%d %s%s: %s -- %s",
            loop,
            "PASS" if new_correct else "FAIL",
            " (FOUND)" if truth_in_candidates else "",
            case["case"],
            confusion_matrix,
        )
        LOG.debug("Results so far: %s", out["totals"])

        if not new_correct:
            LOG.debug(
                "Dot bracket mismatch:\nSEQ: %s\nPRD: %s\nGND: %s",
                case["case"],
                dot_bracket,
                case["truth"],
            )

        out["results"].append(
            {
                # case data
                **case,
                "new_pred": dot_bracket,
                # checks
                "new_correct": new_correct,
                "correct": correct,
                "same_prediction": same_prediction,
                "truth_in_candidates": truth_in_candidates,
                # results
                "results": predictions,
                # scores
                "correct_core_stems": correct_core_stems,
                "confusion_matrix": ", ".join(str(x) for x in confusion_matrix),
            }
        )

    print(json.dumps(out, indent=2))
