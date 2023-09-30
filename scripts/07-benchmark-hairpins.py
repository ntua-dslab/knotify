#!/usr/bin/env python3
#
# Copyright © 2023 Christos Pavlatos, George Rassias, Christos Andrikos,
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
Script:         07-benchmark-hairpins.py
Author:         Angelos Kolaitis <neoaggelos@gmail.com>
Usage:          ./scripts/07-benchmark-hairpins.py cases/hairpins.yaml > out.csv
Description:
    Run benchmark for the cases with hairpins.

    All the results (logfiles and JSON results) are stored in the `results` folder.
    An overview can be seen at `results/out.csv`, which should look like this:

    ```
    ,name,correct,correct_core_stems,truth_in_candidates,count,duration,tp,fp,tn,fn
    0,knotify-base,0,11,0,20,2.498763,292,439,150,338
    1,knotify-hairpins,11,14,13,20,37.933657999999994,628,461,62,68
    2,knotty,3,4,3,20,26.367737,518,366,212,123
    3,ihfold,0,0,0,20,0.683698,484,431,118,186
    4,ihfoldv2,0,4,0,20,1.754078,364,430,110,315
    5,hotknots,5,6,5,20,6.015652,534,418,144,123
    6,ipknot,1,2,1,20,8.782138,540,425,102,152
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
        "ihfold": [
            "--algorithm=ihfold",
        ],
        "ihfoldv2": [
            "--algorithm=ihfoldv2",
        ],
        "hotknots": [
            "--algorithm=hotknots",
        ],
        "ipknot": [
            "--algorithm=ipknot",
        ],
    }

    os.makedirs("results", exist_ok=True)

    for name, args in matrix.items():
        cmd = ["rna_benchmark", "--cases={}".format(cases_yaml), "--verbose", *args]

        print(cmd, file=sys.stderr)
        with open(f"results/{name}.json", "w") as fout, open(
            f"results/{name}.log", "w"
        ) as ferr:
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=ferr)
            stdout = p.stdout.decode()
            print(stdout, file=fout)

        j = json.loads(stdout)
        j["totals"].update(
            {
                "tp": j["totals"]["confusion_matrix"][0],
                "fp": j["totals"]["confusion_matrix"][1],
                "tn": j["totals"]["confusion_matrix"][2],
                "fn": j["totals"]["confusion_matrix"][3],
            }
        )
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
