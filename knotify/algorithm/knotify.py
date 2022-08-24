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

from knotify.algorithm.base import BaseAlgorithm
from knotify.criteria import apply_free_energy_and_stems_criterion
from knotify.energy.base import BaseEnergy
from knotify.energy.vienna import ViennaEnergy
from knotify import hairpin
from knotify.pairalign.base import BasePairAlign
from knotify.parsers.base import BaseParser
from knotify.pairalign.cpairalign import CPairAlign
from knotify.extensions.skip_final_au import SkipFinalAU


class Knotify(BaseAlgorithm):
    def get_results(
        self,
        sequence: str,
        parser: BaseParser,
        skip_final_au: SkipFinalAU = None,
        pairalign: BasePairAlign = None,
        csv: str = None,
        allow_skip_final_au: bool = None,
        max_stem_allow_smaller: int = 1,
        prune_early: bool = False,
        hairpin_grammar: str = None,
        hairpin_allow_ug: bool = False,
        min_hairpin_size: int = hairpin.MIN_HAIRPIN_SIZE,
        min_hairpin_stems: int = hairpin.MIN_HAIRPIN_STEMS,
        max_hairpins_per_loop: int = hairpin.MAX_HAIRPINS_PER_LOOP,
        max_hairpin_bulge: int = hairpin.MAX_HAIRPIN_BULGE,
        energy: BaseEnergy = ViennaEnergy(),
        *args,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Analyze RNA sequence, and predict structure. Return data frame of results
        """
        sequence = sequence.lower()

        pseudoknots = []
        max_size = 0
        for (i, j, left_loop_size, dd_size) in parser(sequence):
            knots = pairalign(sequence, i, j, left_loop_size, dd_size)
            for (dot_bracket, left_loop_stems, right_loop_stems) in knots:
                size = left_loop_stems + right_loop_stems

                if not prune_early or size >= max_size - max_stem_allow_smaller:
                    max_size = max(max_size, size)
                    pseudoknots.append(
                        {
                            "dot_bracket": dot_bracket,
                            "left_loop_stems": left_loop_stems,
                            "right_loop_stems": right_loop_stems,
                            "dd": dd_size,
                        }
                    )

                if not allow_skip_final_au:
                    continue

                knots_without_au = skip_final_au(
                    sequence, dot_bracket, right_loop_stems, left_loop_stems
                )
                if knots_without_au:
                    pseudoknots.extend(
                        {
                            "dot_bracket": d,
                            "right_loop_stems": r,
                            "left_loop_stems": l,
                            "dd": dd_size,
                        }
                        for (d, r, l) in knots_without_au
                    )

        if not pseudoknots:
            pseudoknots = [
                {
                    "dot_bracket": "." * len(sequence),
                    "left_loop_stems": 0,
                    "right_loop_stems": 0,
                    "dd": 0,
                }
            ]

        data = pd.DataFrame(pseudoknots)
        if csv is not None:
            data.to_csv(csv)

        data = apply_free_energy_and_stems_criterion(
            data,
            sequence,
            max_stem_allow_smaller=max_stem_allow_smaller,
            energy=energy,
        )

        if hairpin_grammar is None:
            return data

        data = hairpin.find_hairpins(
            sequence,
            data,
            hairpin_grammar,
            allow_ug=hairpin_allow_ug,
            min_size=min_hairpin_size,
            min_stems=min_hairpin_stems,
            max_bulge=max_hairpin_bulge,
            max_per_loop=max_hairpins_per_loop,
        )
        data = apply_free_energy_and_stems_criterion(
            data,
            sequence,
            max_stem_allow_smaller=max_stem_allow_smaller,
            energy=energy,
        )

        return data
