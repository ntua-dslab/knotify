import json
import yaml

from pseudoknot_detector.main import get_results, argument_parser
from pseudoknot_detector import scoring


def main():
    parser = argument_parser()
    parser.add_argument("--cases", required=True)
    parser.add_argument("--only", type=int, nargs="*")
    parser.add_argument("--correct-stems-slack", type=int, default=0)
    args = parser.parse_args()

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
            "count": len(only),
        },
    }

    for idx in only:
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
        correct = case["pred"] == case["truth"]
        same_prediction = case["pred"] == dot_bracket

        correct_core_stems = scoring.get_correct_core_stems(
            case["truth"],
            dot_bracket,
            slack=args.correct_stems_slack,
        )

        out["totals"]["new_correct"] += new_correct
        out["totals"]["correct"] += correct
        out["totals"]["same_prediction"] += same_prediction
        out["totals"]["correct_core_stems"] += correct_core_stems == 2

        out["results"].append(
            {
                # case data
                **case,
                "new_pred": dot_bracket,
                # checks
                "new_correct": new_correct,
                "correct": correct,
                "same_prediction": same_prediction,
                # results
                "results": predictions,
                # scores
                "correct_core_stems": correct_core_stems,
                "confusion_matrix": ", ".join(
                    str(x)
                    for x in scoring.get_confusion_matrix(
                        case["truth"],
                        dot_bracket,
                        slack=0,
                    )
                ),
            }
        )

    print(json.dumps(out, indent=2))
