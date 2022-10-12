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

from oslo_config import cfg
from knotify.energy.pkenergy import PKEnergy
from knotify.energy.vienna import ViennaEnergy
from knotify.knotify import ENERGY_OPTS

OPTS = [
    cfg.StrOpt("sequence", required=True),
    cfg.StrOpt("dot-bracket", required=True),
]


def main():
    options = cfg.ConfigOpts()
    options.register_cli_opts(OPTS)
    options.register_cli_opts(ENERGY_OPTS)
    options()

    if options.energy == "vienna":
        e = ViennaEnergy()
    elif options.energy == "pkenergy":
        e = PKEnergy(
            options.pkenergy, options.pkenergy_config_dir, options.pkenergy_model
        )

    print(e.eval(options.sequence, options.dot_bracket))


if __name__ == "__main__":
    main()
