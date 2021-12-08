#
# Copyright © 2021 Christos Pavlatos, George Rassias, Christos Andrikos,
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
try:
    from Bio.pairwise2 import align as bioalign
except ModuleNotFoundError:
    from bio.pairwise2 import align as bioalign

# the weights for the pairwise2 globalds alignment algorithm
WEIGHTS = {
    ("a", "u"): 2,
    ("u", "a"): 2,
    ("c", "g"): 2,
    ("g", "c"): 2,
    ("a", "a"): -100,
    ("a", "c"): -100,
    ("a", "g"): -100,
    ("u", "u"): -100,
    ("u", "c"): -100,
    ("u", "g"): 2,
    ("c", "a"): -100,
    ("c", "u"): -100,
    ("c", "c"): -100,
    ("g", "a"): -100,
    ("g", "u"): 2,
    ("g", "g"): -100,
}


def biopairalign(reversed_part, identical_part):
    """Wrapping the bioalign functionality along with the choose function

    :param str to_be_reversed_part: the part that starts with ...
    :param str identical_part: the outer part to be aligned
    :param function choose: the function to apply ordering

    :returns: the most compact alignment, formated based on bio.pairwise policy
    :rtype: list
    """
    # the pairwise alignment score
    best_align_score = 0
    reversed_starts_at = -1
    best_identical_pos = len(identical_part)
    best_pruned_alignment = []
    for i in range(0, len(reversed_part)):

        for j in range(0, len(identical_part)):
            temp = bioalign.globalds(
                reversed_part[i::][::-1],
                identical_part[: len(identical_part) - j :],
                WEIGHTS,
                -3,
                -0.99,
            )

            if len(temp) > 0 and len(temp[0]) > 2 and temp[0][2] >= best_align_score:
                best_align_score = temp[0][2]
                reversed_starts_at = i
                best_identical_pos = j
                best_pruned_alignment = [(el[0][::-1], el[1]) for el in temp]
    identical_ends_at = len(identical_part) - best_identical_pos - 1
    return (reversed_starts_at, identical_ends_at), best_pruned_alignment
