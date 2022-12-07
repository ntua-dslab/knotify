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
from datetime import datetime
import json
import sys
import logging

from oslo_config import cfg
import yaml

from knotify import knotify
from knotify import scoring


LOG = logging.getLogger(__name__)

OPTS = [
    cfg.StrOpt("cases"),
    cfg.ListOpt("only", item_type=cfg.types.Integer()),
    cfg.IntOpt("correct-stems-slack", default=0),
    cfg.BoolOpt("verbose", default=False),
    cfg.BoolOpt("include-candidates", default=False),
    cfg.BoolOpt("include-results", default=True),
]


def main():
    options = knotify.new_options()
    options.register_cli_opts(OPTS)
    options()

    result = run_benchmark(options)

    json.dump(result, sys.stdout, indent=2)


def run_benchmark(options: knotify.ConfigOpts):
    if not options.cases:
        print("Missing required parameter --cases")
        sys.exit(1)

    algorithm, config = knotify.from_options(options)
    if options.verbose:
        logging.basicConfig(level=logging.DEBUG)

        LOG.debug("Command line: %s", sys.argv)
        LOG.debug("Run configuration: %s", config)

    with open(options.cases, "r") as fin:
        cases = yaml.safe_load(fin.read())

    only = options.only or list(range(0, len(cases)))

    out = {
        "results": [],
        "totals": {
            "correct": 0,
            "correct_core_stems": 0,
            "truth_in_candidates": 0,
            "count": len(only),
            "duration": 0,
            "confusion_matrix": [0, 0, 0, 0],
        },
    }

    for loop, idx in enumerate(only):
        case = cases[idx]
        start = datetime.now()
        results = algorithm.get_results(case["case"].lower(), **config)
        duration = datetime.now() - start
        candidates = results[["dot_bracket", "energy", "stems"]].to_dict(
            orient="records"
        )
        dot_bracket = candidates[0]["dot_bracket"]
        energy = candidates[0]["energy"]

        correct = case["truth"] == dot_bracket

        try:
            correct_core_stems = scoring.get_correct_core_stems(
                case["truth"], dot_bracket, slack=options.correct_stems_slack
            )
        except ValueError as e:
            LOG.exception("Failed to retrieve number of correct core stems: %s", e)
            correct_core_stems = 0

        confusion_matrix = scoring.get_confusion_matrix(
            case["truth"], dot_bracket, slack=0
        )
        truth_in_candidates = case["truth"] in [p["dot_bracket"] for p in candidates]

        out["totals"]["correct"] += correct
        out["totals"]["correct_core_stems"] += correct_core_stems == 2
        out["totals"]["truth_in_candidates"] += truth_in_candidates
        out["totals"]["duration"] += duration.total_seconds()

        out["totals"]["confusion_matrix"][0] += confusion_matrix[0]
        out["totals"]["confusion_matrix"][1] += confusion_matrix[1]
        out["totals"]["confusion_matrix"][2] += confusion_matrix[2]
        out["totals"]["confusion_matrix"][3] += confusion_matrix[3]

        LOG.info(
            "%d %s%s: %s -- %.2f seconds -- %s",
            loop,
            "PASS" if correct else "FAIL",
            " (FOUND)" if truth_in_candidates else "",
            case["case"],
            duration.total_seconds(),
            confusion_matrix,
        )
        LOG.debug("Results so far: %s", out["totals"])

        if not correct:
            LOG.debug(
                "Dot bracket mismatch:\nSEQ: %s\nPRD: %s\nGND: %s",
                case["case"],
                dot_bracket,
                case["truth"],
            )

        item = {
            # case data
            **case,
            "pred": dot_bracket,
            # checks
            "correct": correct,
            "truth_in_candidates": truth_in_candidates,
            # scores
            "correct_core_stems": correct_core_stems,
            "confusion_matrix": ", ".join(str(x) for x in confusion_matrix),
            "duration": duration.total_seconds(),
            "energy": energy,
        }

        if options.include_candidates:
            item["candidates"] = candidates

        if options.include_results:
            out["results"].append(item)

    return out
