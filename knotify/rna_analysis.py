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
from knotify.models import Pknot
from knotify.pairing import (
    get_left_stem_aligned_indices,
    get_right_stem_aligned_indices,
)
from knotify.parsers.base import BaseParser


def replace_str_at_index(string, index, char):
    return string[:index] + char + string[index + 1 :]


class StringAnalyser(object):
    def __init__(self, input_string: str, parser: BaseParser = None):
        """Basic rna analysis class

        :param str input: the input string to be analyzed
        :param BaseParser parser: RNA parser instance
        """
        self.sequence = input_string
        self.parser = parser

    def get_pseudoknots_without_au(self, sequence: str, knot: Pknot):
        """
        return list of pseudoknots without AU as last loop stem.
        """
        # TODO(akolaitis): this is quite dirty. clean up code
        # and remove duplication
        l = []
        removed = 0
        if all(knot.right_loop_stems):
            inner_idx = knot.right_loop_stems[0][0]
            outer_idx = knot.right_loop_stems[1][-1]
            if set([sequence[inner_idx], sequence[outer_idx]]) == set("au"):
                removed += 1
                knot2 = knot.to_dict()
                knot2["right_loop_stems"] -= 1
                knot2["dot_bracket"] = replace_str_at_index(
                    knot2["dot_bracket"], inner_idx, "."
                )
                knot2["dot_bracket"] = replace_str_at_index(
                    knot2["dot_bracket"], outer_idx, "."
                )
                l.append(knot2)

        if all(knot.left_loop_stems):
            inner_idx = knot.left_loop_stems[0][-1]
            outer_idx = knot.left_loop_stems[1][0]
            if set([sequence[inner_idx], sequence[outer_idx]]) == set("au"):
                removed += 1
                knot2 = knot.to_dict()
                knot2["left_loop_stems"] -= 1
                knot2["dot_bracket"] = replace_str_at_index(
                    knot2["dot_bracket"], inner_idx, "."
                )
                knot2["dot_bracket"] = replace_str_at_index(
                    knot2["dot_bracket"], outer_idx, "."
                )

                l.append(knot2)

        if removed == 2:
            knot2 = knot.to_dict()
            knot2["left_loop_stems"] -= 1
            knot2["right_loop_stems"] -= 1

            for idx in [
                knot.right_loop_stems[0][0],
                knot.right_loop_stems[1][-1],
                knot.left_loop_stems[0][-1],
                knot.left_loop_stems[1][0],
            ]:
                knot2["dot_bracket"] = replace_str_at_index(
                    knot2["dot_bracket"], idx, "."
                )

            l.append(knot2)

        return l

    def get_pseudoknots(
        self, max_stem_allow_smaller=2, prune_early=False, allow_skip_final_au=False
    ):
        pseudoknots = []
        max_size = 0
        for line in self.parser.detect_pseudoknots(self.sequence):
            i, j, left_loop_size, dd_size = line
            knot = self._get_pseudoknots_for_tree(i, j, left_loop_size, dd_size)
            size = len(knot.right_loop_stems[0]) + len(knot.left_loop_stems[0])
            if knot and (not prune_early or size >= max_size - max_stem_allow_smaller):
                max_size = max(max_size, size)
                pseudoknots.append(knot.to_dict())

                if not allow_skip_final_au:
                    continue

                knots_without_au = self.get_pseudoknots_without_au(self.sequence, knot)
                if knots_without_au:
                    pseudoknots.extend(knots_without_au)

        return pseudoknots

    def _get_pseudoknots_for_tree(self, i, j, left_loop_size, dd_size):
        """Returns the pseudoknot for any given substring

        :param int i: the start of the sliding window
        :param int j: the end of the sliding window
        :param yaep_serialized tree: the tree to look pseudoknots on, generated
                by c code

        :returns: a pseudoknot
        :rtype: Pknot
        """
        pknot = Pknot(len=len(self.sequence))
        pknot.left_core_indices = (i, i + left_loop_size + dd_size + 2)

        pknot.right_core_indices = (i + left_loop_size + 1, i + j - 1)
        # ========================================================================
        # get strings to align with pairwise
        _, left_stem_indices = get_left_stem_aligned_indices(
            self.sequence,
            pknot.get_left_inner_potential(),
            pknot.get_left_outer_potential(),
        )

        pknot.left_loop_stems = left_stem_indices

        _, right_stem_indices = get_right_stem_aligned_indices(
            self.sequence,
            pknot.get_right_inner_potential(),
            pknot.get_right_outer_potential(),
        )

        pknot.right_loop_stems = right_stem_indices

        return pknot
