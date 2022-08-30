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


@pytest.mark.parametrize("parser", ["bruteforce", "yaep"])
@pytest.mark.parametrize(
    "name, sequence, candidate, overrides",
    [
        (
            "without final au",
            "GGGAAACGAGCCAAGUGGCGCCGACCACUUAAAAACACCGGAA",
            ".............(((((..[[[.))))).........]]]..",
            {"allow_skip_final_au": True},
        ),
        *[
            (
                "without final au (both sides)",
                "GGGAAACGAGCCAAGUGGCUCCGACCACUUAAAAACACCGGAA",
                candidate,
                {
                    "allow_skip_final_au": True,
                    "max_stem_allow_smaller": 3,
                },
            )
            for candidate in [
                ".............(((((..[[[.))))).........]]]..",  # drop both
                ".............(((((.[[[[.))))).........]]]].",  # drop left
                "............((((((..[[[.))))))........]]]..",  # drop right
                "............((((((.[[[[.))))))........]]]].",  # include all
            ]
        ],
    ],
)
def test_exploration(parser, name, sequence, candidate, overrides):
    opts = knotify.new_options()

    opts.parser = parser
    for key, value in overrides.items():
        opts.set_override(key, value)

    algorithm, config = knotify.from_options(opts)

    assert candidate in algorithm.get_results(sequence, **config).dot_bracket.unique()
