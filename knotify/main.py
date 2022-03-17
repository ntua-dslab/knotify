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
import sys

from oslo_config import cfg

from knotify import knotify

OPTS = [
    cfg.StrOpt("sequence"),
]


def main():
    options = knotify.new_options()
    options.register_cli_opts(OPTS)
    options()

    if not options.sequence:
        print("Missing required parameter --sequence")
        sys.exit(1)

    algorithm, config = knotify.from_options(options)

    start = datetime.now()
    results = algorithm.get_results(sequence=options.sequence.lower(), **config)
    duration = datetime.now() - start

    if options.results_csv:
        results[["dot_bracket", "stems", "energy"]].to_csv(
            options.results_csv, index=None
        )

    chosen = results.loc[0]

    print("Sequence: ", options.sequence)
    print("Structure:", chosen.dot_bracket)
    print("Energy:", chosen.energy)
    print("Duration:", duration.total_seconds(), "s")


if __name__ == "__main__":
    main()
