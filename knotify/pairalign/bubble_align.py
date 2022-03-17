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


def string_align(str1, str2):
    """Given two strings it aligns them in
    a custom way to minimize holes and
    hole lengths
    """

    # length align the two strings
    common_length()

    align1 = []
    align2 = []
    for i, j, k in enumerate(zip(str1, str2)):
        if i == j:
            align1.append(i)
            align2.append(j)
        else:
            next_candidates = [
                (str1[k:], str2[k + 1 :]),
                (str1[k + 1 :], str2[k:]),
                (str1[k + 1 :], str2[k + 1 :]),
            ]
            partial_alignments = [
                string_align(*substrings) for substrings in next_candidates
            ]

            return (
                "".join(align1) + partial_alignments[0][1],
                "".join(align2) + partial_alignments[0][1],
            )


# split the string into segments
