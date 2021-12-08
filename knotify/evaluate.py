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
import json
import yaml

from knotify.main import get_results, argument_parser
from knotify import scoring


def main():
    parser = argument_parser()
    parser.add_argument("--cases", required=True)
    parser.add_argument("--correct-stems-slack", type=int, default=0)
    args = parser.parse_args()

    with open(args.cases, "r") as fin:
        cases = yaml.safe_load(fin.read())

    results, correct, no_pseudoknot = [], 0, 0
    for case in cases:
        try:
            correct_core_stems = scoring.get_correct_core_stems(
                case["truth"],
                case["pred"],
                slack=args.correct_stems_slack,
            )
        except (ValueError, IndexError):
            correct_core_stems = "invalid"
            no_pseudoknot += 1

        ok = scoring.find_matches(case["truth"]) == scoring.find_matches(case["pred"])
        if ok:
            correct += 1

        results.append(
            {
                # case data
                **case,
                "correct_core_stems": correct_core_stems,
                "confusion_matrix": ", ".join(
                    str(x)
                    for x in scoring.get_confusion_matrix(case["truth"], case["pred"])
                ),
                "correct": ok,
            }
        )

    print(
        yaml.dump(
            {
                "config": {
                    "correct_core_stems_slack": args.correct_stems_slack,
                },
                "results": results,
                "totals": {
                    "total": len(cases),
                    "correct": correct,
                    "no_pseudoknot": no_pseudoknot,
                },
            }
        )
    )


if __name__ == "__main__":
    main()
