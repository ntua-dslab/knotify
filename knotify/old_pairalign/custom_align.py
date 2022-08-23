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
from operator import itemgetter

# this should be a global settings TODO: later on
COMBINATIONS = ["au", "ua", "gc", "cg", "ug", "gu"]


def align_consecutive(str1, str2):
    no_iter = min(len(str1), len(str2))
    str1_aligned = []
    str2_aligned = []
    i = 0

    # get the first part
    while i < no_iter and (str1[i] + str2[i]) in COMBINATIONS:
        str1_aligned.append(str1[i])
        str2_aligned.append(str2[i])
        i += 1
    # hint i is the open string limit, not an index
    return i, "".join(str1_aligned), "".join(str2_aligned)


def bubble_align(str1, str2):
    """Returns the OPEN ends for str1, str2 along with the bubbled alignments
    OPEN end is the first not matching character position
    """

    open_end, str1_aligned_prefix, str2_aligned_prefix = align_consecutive(str1, str2)

    # if no potential gaps...
    if open_end >= len(str1) or open_end >= len(str2):
        return open_end, open_end, str1_aligned_prefix, str2_aligned_prefix
    elif open_end < 1:
        return 0, 0, "", ""

    next_candidates = [
        (str1[open_end:], str2[open_end + 1 :]),
        (str1[open_end + 1 :], str2[open_end:]),
        (str1[open_end + 1 :], str2[open_end + 1 :]),
    ]

    candidate_results = [
        align_consecutive(*next_candidate) for next_candidate in next_candidates
    ]

    relative_open_end, str1_aligned_suffix, str2_aligned_suffix = max(
        candidate_results, key=itemgetter(0)
    )

    if relative_open_end == 0:
        return open_end, open_end, str1_aligned_prefix, str2_aligned_prefix

    # check for the suffix rule
    rule_index = candidate_results.index(max(candidate_results, key=itemgetter(0)))

    # now we have distinct open end for the two strings
    str1_open_end = open_end + relative_open_end
    str2_open_end = open_end + relative_open_end

    # based on the rule index we should introduce the corresponding bubble
    if rule_index == 0:
        str1_aligned = str1_aligned_prefix + str1_aligned_suffix
        str2_aligned = str2_aligned_prefix + "-" + str2_aligned_suffix
    elif rule_index == 1:
        str1_aligned = str1_aligned_prefix + "-" + str1_aligned_suffix
        str2_aligned = str2_aligned_prefix + str2_aligned_suffix
    else:
        str1_aligned = str1_aligned_prefix + "-" + str1_aligned_suffix
        str2_aligned = str2_aligned_prefix + "-" + str2_aligned_suffix

    # str1_open_end, str2_open_end are open end limits not indices
    return str1_open_end, str2_open_end, str1_aligned, str2_aligned


def custom_pairalign(reversed_part, identical_part):

    """Reversed part is the one to be reversed instead of already being
    reversed
    """
    (
        reversed_ends_at,
        identical_ends_at,
        reversed_aligned,
        identical_aligned,
    ) = bubble_align(reversed_part[::-1], identical_part)
    max_align_len = max(len(reversed_aligned), len(identical_aligned))
    min_align_len = min(len(reversed_aligned), len(identical_aligned))

    # add any missing bubble in the beggining or the end of each string
    # to fully align them

    # the return values should be INDICES of the first and last elment respectiverly
    #  bellow
    return (len(reversed_part) - len(reversed_aligned), len(identical_aligned) - 1), (
        "".join(reversed_aligned[::-1]),
        "".join(identical_aligned),
    )


# for testing:
# custom_pairalign("uaggggugucagggcucug", "cgagccccggcuaacucuaucugua")
#     ...:
# Out[55]: (12, 6, 'gggcu-ug', 'cgagccc')
