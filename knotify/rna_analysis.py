from time import time
import multiprocessing

from knotify.models import Pknot, OUTER, INNER
from knotify.pairing import (
    get_left_stem_aligned_indices,
    get_right_stem_aligned_indices,
)
from knotify.rna_parser import PseudoknotDetector

MAX_LOOP_SIZE = 100  # max number of unpaired bases belonging to a loop
MIN_LOOP_SIZE = 1  # min number of unpaired bases belonging to a loop


def replace_str_at_index(string, index, char):
    return string[:index] + char + string[index + 1 :]


class StringAnalyser(object):
    def __init__(
        self,
        input_string,
        grammar=None,
        max_loop_size=MAX_LOOP_SIZE,
        max_dd_size=2,
        min_dd_size=0,
        allow_ug=True,
    ):
        """Basic rna analysis class

        :param str input: the input string to be analyzed
        :param int max_loop_size: the maximum length of any of the two loops of
            the pseudoknot
        :param str grammar: the analysis grammar
        """
        self._string = input_string
        self._max_window_size = 2 * MAX_LOOP_SIZE + 4
        self._min_window_size = 2 * MIN_LOOP_SIZE + 4
        self.parser = PseudoknotDetector(
            grammar=grammar,
            max_dd_size=max_dd_size,
            allow_ug=allow_ug,
            min_dd_size=min_dd_size,
        )
        self.window_boundaries = []

    def get_window_boundaries(self):
        """Generates all window boundaries and store them to the class

        :return: the StringAnalyser object (pipeline)

        :rtype: object
        """
        # TODO (akolaitis): most of the heavy lifting here could be easily
        # handed over to the c code. this means that the C code should also
        # return the left and right loop index (apart from the left loop
        # boundary and the dd size). currently, these are are implied by the
        # "i" and "j" indices of the tree.
        #
        # A naive solution would be to include the "i" and "j" in each result
        # but that greatly increases the required memory for the intermediate
        # results, e.g.
        # (for i1, j1)
        # "5, 7" --> "i1, j1, 5, 7"
        # "5, 8" --> "i1, j1, 5, 8"
        # (for i2, j2)
        # "4, 6" --> "i2, j2, 4, 6"
        # "4, 7" --> "i2, j2, 4, 7"
        #
        # Consider an improved version where results are grouped by "i" and "j",
        # e.g. with the following output structure (different separators):
        # "i1,j1:5,7|5,8"
        # "i2,j2:4,7|4,7"
        window_boundaries = []
        for i in range(0, len(self._string) - self._min_window_size):
            window_size = self._max_window_size
            if (i + self._max_window_size) > len(self._string):
                window_size = len(self._string) - i
            for j in range(window_size, self._min_window_size, -1):
                window_boundaries.append((i, j))
        self.window_boundaries = window_boundaries

        return self

    def generate_trees_in_parallel(self):
        start = time()
        with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
            trees = p.map(
                self.parser.detect_pseudoknots,
                [self._string[s : s + l] for s, l in self.window_boundaries],
            )

        self.trees = trees
        return self

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
        self,
        max_stem_allow_smaller=2,
        prune_early=False,
        allow_skip_final_au=False,
    ):
        # TODO(akolaitis) 1: move logic out of the for loop in separate function
        # and then run this loop within with detect_pseudoknots

        # TODO(akolaitis) 2: this can be done in pandas as well
        pseudoknots = []
        max_size = 0
        for ((i, j), trees) in zip(self.window_boundaries, self.trees):
            for tree in trees:
                if tree:
                    knot = self._get_pseudoknots_for_tree(i, j, tree)
                    size = len(knot.right_loop_stems[0]) + len(knot.left_loop_stems[0])
                    if knot and (
                        not prune_early or size >= max_size - max_stem_allow_smaller
                    ):
                        max_size = max(max_size, size)
                        pseudoknots.append(knot.to_dict())

                        if not allow_skip_final_au:
                            continue

                        knots_without_au = self.get_pseudoknots_without_au(
                            self._string, knot
                        )
                        if knots_without_au:
                            pseudoknots.extend(knots_without_au)

        return pseudoknots

    def left_loop_boundary(self, tree):
        return int(tree.split(",")[0])

    def dd(self, tree):
        return int(tree.split(",")[1])

    def _get_pseudoknots_for_tree(self, i, j, tree):
        """Returns the pseudoknot for any given substring

        :param int i: the start of the sliding window
        :param int j: the end of the sliding window
        :param yaep_serialized tree: the tree to look pseudoknots on, generated
                by c code

        :returns: a pseudoknot
        :rtype: Pknot
        """
        pknot = Pknot(len=len(self._string))
        pknot.left_core_indices = (
            i,
            i + self.left_loop_boundary(tree) + self.dd(tree) + 2,
        )
        pknot.i_j = (i, j, self.left_loop_boundary(tree))

        pknot.right_core_indices = (i + self.left_loop_boundary(tree) + 1, i + j - 1)

        pknot.tree = tree
        pknot.ij = "{i}, {j}".format(i=i, j=j)
        # ========================================================================
        # get strings to align with pairwise
        _, left_stem_indices = get_left_stem_aligned_indices(
            self._string,
            pknot.get_left_inner_potential(),
            pknot.get_left_outer_potential(),
        )

        pknot.left_loop_stems = left_stem_indices

        _, right_stem_indices = get_right_stem_aligned_indices(
            self._string,
            pknot.get_right_inner_potential(),
            pknot.get_right_outer_potential(),
        )

        pknot.right_loop_stems = right_stem_indices

        # check if pseudoknot is out of bounds or not
        # if left_loop_size > MAX_LOOP_SIZE or right_loop_size > MAX_LOOP_SIZE:
        #     return None

        return pknot

    def _get_no_stems(self, str1, str2):
        """Given two strings, it returns the number of consequent stems
        counting from the beggining of the two

        :param str str1: the first string
        :param str str2: the second string

        :return: the number of stems
        :rtype: int
        """
        count = 0
        no_iter = min(len(str1), len(str2))
        for i in range(0, no_iter):
            if (str1[i] + str2[i]) in ["au", "ua", "gc", "cg", "ug", "gu"]:
                count += 1
            else:
                break
        return count


def plot_trees_distribution_to_file(
    starting_point_or_length, stem_sums, file_name=None
):
    file_name = "trees_scattered_distribution"
    xs = starting_point_or_length
    ys = stem_sums
    from matplotlib import pyplot as plt

    plt.scatter(xs, ys)
    plt.savefig("file_name.{}".format("png"))
    return 0


# ====================== helper functions ===============================


def analyze_that(sequence):
    """The overall exposed functionality"""
    sequence = sequence.lower()
    analyzer = StringAnalyser(sequence)
    start_timestamp = time()

    analyzer.get_window_boundaries().generate_trees_in_parallel()

    # hack to fix the no ambiguous or term trees issue
    # TODO: the following line may be omitted
    analyzer._pos_trees = [tree for tree in analyzer._pos_trees if tree[2] != ""]

    pseudoknots = analyzer.get_pseudoknots()
    print("elapsed time: {} seconds".format(time() - start_timestamp))
    stem_sums = [
        len(pseudoknot.left_loop_stems[0]) + len(pseudoknot.right_loop_stems[0])
        for pseudoknot in pseudoknots
    ]

    return stem_sums, pseudoknots
    min_en_indices = [i for (i, j) in enumerate(stem_sums) if j == max(stem_sums)]

    for min_en_index in min_en_indices:
        print(analyzer._pos_trees[min_en_index])
        dot = DotBracket(**pseudoknots[min_en_index], original_string=sequence).print()

        k = DeflatedPseudoknot(dot, **pseudoknots[min_en_index])
        output = "{info}\n{input_str}\n{dot}".format(
            info=pseudoknots[min_en_index], input_str=sequence, dot=dot
        )

        print(output)
        end_timestamp = time()
        print(
            "\nThe output has been calculated in: {} seconds".format(
                end_timestamp - start_timestamp
            )
        )


if __name__ == "__main__":
    # input1 = "UGUCGGUGGUCAUUGCGGAGGGGGAACGCCCGGUCCCAUCCCGAACCCGGAAGCUAAGCCCUCCAGCGCCGAUGGUACUGCACUCGCCAGGGUGUGGGAGAGUAGGUCGCCGCCGACA".lower()
    input1 = "cgugaaggcuacgauagugcc".lower()
    # input1 = "AGCUGCCCUUGGGUUUUACUCCUUGAACCCUUCGGAAGAACUCUUUGGAGUUCGUACCAGUACCUCACAUAGUGAGGUAAUAAGACUGGUGGGCAGCGCCUAGUCGAAAGACUAGGUGAUCUCUAAGGAGACCA".lower()

    analyze_that(input1)
