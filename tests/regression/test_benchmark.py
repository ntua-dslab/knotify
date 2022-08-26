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
import pytest

from knotify import knotify
from knotify import benchmark


@pytest.mark.parametrize(
    "name,options,expected",
    [
        (
            "2022-basic",
            {
                "cases": "./cases/cases.yaml",
                "max_dd_size": 2,
                "max_stem_allow_smaller": 2,
                "allow_ug": True,
                "prune_early": True,
                "parser": "bruteforce",
                "energy": "vienna",
                "allow_skip_final_au": True,
            },
            {
                "correct": 11,
                "correct_core_stems": 15,
                "truth_in_candidates": 14,
            },
        ),
        (
            "2022-aiai",
            {
                "cases": "./cases/new.yaml",
                "max_dd_size": 2,
                "max_stem_allow_smaller": 1,
                "allow_ug": True,
                "prune_early": True,
                "parser": "bruteforce",
                "energy": "pkenergy",
                "min_window_size_ratio": 0.26,
                "max_window_size_ratio": 0.87,
                "min_window_size": 0,
                "max_window_size": 0,
            },
            {
                "correct": 81,
                "correct_core_stems": 141,
                "truth_in_candidates": 106,
            },
        ),
    ],
)
def test_benchmark(name, options, expected):
    opts = knotify.new_options()
    opts.register_cli_opts(benchmark.OPTS)

    for key, value in options.items():
        opts.set_override(key, value)

    result = benchmark.run_benchmark(opts)

    for key, value in expected.items():
        assert result["totals"][key] == value
