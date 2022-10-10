#!/bin/bash

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

# Script:         05-benchmark-bulges-results.sh
# Author:         Angelos Kolaitis <neoaggelos@gmail.com>
# Usage:          ./scripts/05-benchmark-bulges-results.sh cases/new.yaml
# Description:
#     Run benchmark tests for IEEE 2022.
#     Compare:
#     - knotify-old: Old knotify implementation (no support for bulges)
#     - knotify-bulges: Knotify implementation with support for bulges
#     - knotty: Knotty implementation, with support for bulges (with baseline configuration)

rna_benchmark --cases cases/new.yaml --max-dd-size=2 --parser=bruteforce --allow-ug --max-stem-allow-smaller=2 --energy=pkenergy --allow-skip-final-au --pairalign=consecutive --verbose > knotify-old.json 2> knotify-old.log

rna_benchmark --cases cases/new.yaml --max-dd-size=2 --parser=bruteforce --allow-ug --max-stem-allow-smaller=2 --energy=pkenergy --allow-skip-final-au --pairalign=bulges --max-bulge-size=3 --min-stems-after-bulge=2 --verbose > knotify-bulges.json 2> knotify-bulges.log

rna_benchmark --cases cases/new.yaml --algorithm knotty --verbose > knotty.json 2> knotty.log
