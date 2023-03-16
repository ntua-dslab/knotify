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
Script:         07-benchmark-bulges-initial.py
Author:         Angelos Kolaitis <neoaggelos@gmail.com>
Usage:          ./scripts/04-benchmark-bulges-initial.py cases/hairpins.yaml > out.csv
Description:
    Run benchmark for the cases with bulges.

    The output is printed as a CSV file, like so:

    ```
    ,name,correct,correct_core_stems,truth_in_candidates,count,duration
    0,knotty,4,15,4,35,36.922723999999995
    1,knotify-old,0,4,0,35,26.601583
    2,knotify-new,3,13,4,35,18.831148000000002
    ```
"""

import argparse
import json
import subprocess
import sys

import os

import pandas as pd


def run_benchmarks(cases_yaml: str):
    records = []

    matrix = {
        "knotify-base": [
            "--max-dd-size=2",
            "--parser=bruteforce",
            "--allow-ug",
            "--max-stem-allow-smaller=2",
            "--energy=pkenergy",
            "--pairalign=consecutive",
        ],
        "knotify-hairpins": [
            "--max-dd-size=2",
            "--parser=bruteforce",
            "--allow-ug",
            "--max-stem-allow-smaller=2",
            "--energy=pkenergy",
            "--pairalign=consecutive",
            "--hairpin-grammar=./libhairpin.so",
        ],
        "knotty": [
            "--algorithm=knotty",
        ],
    }

    os.makedirs("results", exist_ok=True)

    for name, args in matrix.items():
        cmd = ["rna_benchmark", "--cases={}".format(cases_yaml), "--verbose", *args]

        print(cmd, file=sys.stderr)
        with open(f"results/{name}.json", "w") as fout, open(f"results/{name}.log", "w") as ferr:
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=ferr)
            stdout = p.stdout.decode()
            print(stdout, file=fout)

        j = json.loads(stdout)
        j["totals"].update({
            "tp": j["totals"]["confusion_matrix"][0],
            "fp": j["totals"]["confusion_matrix"][1],
            "tn": j["totals"]["confusion_matrix"][2],
            "fn": j["totals"]["confusion_matrix"][3],
        })
        del j["totals"]["confusion_matrix"]

        records.append({"name": name, **j["totals"]})

    with open("results/out.csv", "w") as fout:
        pd.DataFrame(records).to_csv(fout)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cases")

    args = parser.parse_args()
    return run_benchmarks(args.cases)


if __name__ == "__main__":
    main()
