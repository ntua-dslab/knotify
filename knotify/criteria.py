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
import pandas as pd

from knotify.energy.base import BaseEnergy


def apply_free_energy_and_stems_criterion(
    data: pd.DataFrame,
    sequence: str,
    max_stem_allow_smaller: int,
    energy: BaseEnergy,
):
    """
    Returns the best result for a Pandas data frame.

    Expected shape for the Pandas DataFrame rows:
    {
        "left_loop_stems": len(knot.left_loop_stems[0]),
        "right_loop_stems": len(knot.right_loop_stems[0]),
        "dot_bracket": visualize_knot(knot),
        "dd": str(len(knot.get_dd_seg()))
    }
    """
    # TODO(akolaitis): Consider supporting "strategies" where smaller number of stems
    # are allowed but are discarded based on energy. Also consider favoring pseudoknots
    # with larger loop sizes, which allows for more hairpin/bulge options.
    #
    # Note that calculating energy for all dot brackets is a very expensive
    # operation. However, we have seen dot brackets with few stems that give
    # very good energy results (and are part of the ground truth, because of
    # bulges and/or hairpins in the pseudoknot loops). For an example, see case
    # with identifier "knotify/case/2021/paper1/bulges/1".
    #
    # As a future reference, benchmarking energy calculation:
    # Sequence: guuucuccuaucgccaucuggaugggauuuaagagacuuaugccaaacucuuugaguuugaugccaauucaguuuucgucugaauugaugcccgaaaggauggc    #  noqa
    # calculating energy for 85554 rows, finished after 7.016129016876221 sec
    # calculating energy for 2949 rows, finished after 0.2533254623413086 sec

    # max stems
    data["stems"] = data["left_loop_stems"] + data["right_loop_stems"]
    data["real_stems"] = data["dot_bracket"].apply(lambda r: sum(x != "." for x in r))
    data = data[data["stems"] >= data["stems"].max() - max_stem_allow_smaller].reset_index()

    # min energy
    data["energy"] = data["dot_bracket"].apply(lambda r: energy.eval(sequence, r))
    data.sort_values(["energy", "real_stems", "dd"], ascending=(True, False, True), inplace=True)

    data = data.reset_index()

    return data
