#!/usr/bin/env python3

"""
Script:         00-example.py
Author:         Angelos Kolaitis <neoaggelos@gmail.com>
Usage:          ./00-example.py
Description:
    This script gives an example of using Knotify directly from Python code.
"""

from datetime import datetime

# import knotify
from knotify import knotify


def main():
    # initialize knotify options
    options = knotify.new_options()

    # override any options. you will most likely need to set the C library paths.
    # refer to knotify.*_OPTS for defaults.
    #
    # For example,
    #
    # options.yaep_library_path = "/opt/knotify/lib/libpseudoknot.so"

    # initialize the knotify engine from options
    algorithm, config = knotify.from_options(options)

    # run knotify and get list of results
    start = datetime.now()
    results = algorithm.get_results("CGAAUCUCAAGCAAUCAAGCAUUCUACUUCUAUUGCA", **config)
    duration = datetime.now() - start

    # first result is our choice, the rest are candidates
    result = results.loc[0]

    print("Dot Bracket:", result.dot_bracket)
    print("Energy:", result.energy)
    print("Duration:", duration.total_seconds())


if __name__ == "__main__":
    main()
