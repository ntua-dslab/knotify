#!/usr/bin/env python3

"""
Script:         01-pseudoknot-stats.py
Author:         Angelos Kolaitis <neoaggelos@gmail.com>
Usage:          ./01-pseudoknot-stats.py <cases.yaml>
Description:
    Given a file containing RNA sequences and their respective ground truths,
    calculate the following metrics, and print them in csv format to stdout:

    | Metric    | Description |
    | --------- | ----------- |
    | num_stems | Number of stems in the ground truth |
    | pknot_len | Pseudoknot length (right outer core - left outer core) |
    | unused_left | Number of unused bases (up to left outer core) |
    | unused_right | Number of unused bases (up to right outer core) |

    For example, consider the following input:

    ```yaml
    # cases.yaml
    - case:  CGAAUCUCAAGCAAUCAAGCAUUCUACUUCUAUUGCA
      truth: .((((((...[[[[[..)).)))).......]]]]].
    ```

    Run the script with:
    ```bash
    ./scripts/01-pseudoknot-stats.py ./cases.yaml > output.csv
    ```

    The output for this particular case would be:

    ```csv
    # output.csv
    ,truth,case_size,num_stems,pknot_size,left_unused,right_unused
    0,.((((((...[[[[[..)).)))).......]]]]].,37,22,25,6,7
    ```

    The detailed explanation follows:

    ```
    0         1         2         3
    0123456789012345678901234567890123456 | length = 37
    CGAAUCUCAAGCAAUCAAGCAUUCUACUUCUAUUGCA
    .((((((...[[[[[..)).)))).......]]]]].
          ^                        ^
    ......(.......[..).............].....
    ```

    Left outer core is at index 6, and right outer core is at index 31.

    num_stems = 6 + 5 + 6 + 5 = 22
    pknot_size = 31 - 6 = 25
    left_unused = 6
    right_unused = 37 - 31 - 1 = 5
"""

import argparse
import sys

import yaml
import pandas as pd


def print_stats(cases_yaml: str):
    with open(cases_yaml) as fin:
        cases = yaml.safe_load(fin)

    results = []
    for truth in [c["truth"] for c in cases]:
        case_length = len(truth)
        first_right_stem = truth.find("[")
        left = truth[:first_right_stem].rfind("(")
        right = truth.find("]")
        results.append(
            {
                "truth": truth,
                "case_size": case_length,
                "num_stems": case_length - sum(t == "." for t in truth),
                "pknot_size": right - left,
                "left_unused": left,
                "right_unused": case_length - right - 1,
            }
        )

    pd.DataFrame(results).to_csv(sys.stdout)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cases")

    args = parser.parse_args()
    return print_stats(args.cases)


if __name__ == "__main__":
    main()
