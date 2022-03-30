#!/usr/bin/env python3
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

"""
Script:         03-benchmark-window-sizes-b.py
Author:         Angelos Kolaitis <neoaggelos@gmail.com>
Usage:          ./scripts/03-benchmark-window-sizes-b.py cases/new.yaml > out.csv
Description:
    Similar output with script 02-benchmark-window-sizes.py.
    This is adding a few more test cases with finer-grained tuning on the window sizes.
"""

import argparse
import json
import subprocess
import sys

import pandas as pd


def run_benchmarks(cases_yaml: str):
    records = []

    for parser in ["bruteforce", "yaep"]:
        for minw, maxw in [
            (0.31, 0.82),
            (0.26, 0.87),
            (0.21, 0.92),
            (0.16, 0.97),
            (0.06, 1),
        ]:
            cmd = [
                "rna_benchmark",
                "--cases={}".format(cases_yaml),
                "--max-dd-size=2",
                "--max-stem-allow-smaller=1",
                "--allow-ug",
                "--prune-early",
                "--parser={}".format(parser),
                "--energy=pkenergy",
                "--min-window-size-ratio={}".format(minw),
                "--max-window-size-ratio={}".format(maxw),
                "--min-window-size=0",
                "--max-window-size=0",
                "--noinclude-results",
            ]

            print(cmd, file=sys.stderr)
            stdout = subprocess.check_output(cmd)
            j = json.loads(stdout)

            records.append({"minw": minw, "maxw": maxw, "parser": parser, **j["totals"]})

    pd.DataFrame(records).to_csv(sys.stdout)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cases")

    args = parser.parse_args()
    return run_benchmarks(args.cases)


if __name__ == "__main__":
    main()
