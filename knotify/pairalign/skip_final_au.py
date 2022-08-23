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
import ctypes


class SkipFinalAU:
    """
    Accepts a dot bracket that has been decorated with the loop stems. If the
    final stem of the right or the left loop is an AU pair, return an extra dot
    bracket without that final stem.

    The implementation is done in C code in pairalign/skipfinalau.c

    Example:

    AAGCCUUG   --- sequence
    ((([)))]   --- dot bracket

    Will return the following dot brackets:

    ((([)))]   --- original dot bracket
    .(([)).]   --- dot bracket without last AU pair on left loop
    """

    def __init__(self, library_path: str, *args, **kwargs):
        super(SkipFinalAU, self).__init__(*args, **kwargs)

        self.lib = ctypes.CDLL(library_path)

    def pairalign(
        self,
        sequence: str,
        dot_bracket: str,
        left_loop_stems: int,
        right_loop_stems: int,
    ) -> list:
        """
        Accept a dot bracket that has been decorated with loop stems. Return a list
        of (dot_bracket, right_loop_stems, left_loop_stems) tuples if it is possible
        to trim the last AU pair from the loop stems.
        """
        results = []

        def add_result(dot_bracket, left_loop_stems, right_loop_stems):
            results.append(
                (
                    dot_bracket.decode(),
                    left_loop_stems,
                    right_loop_stems,
                )
            )

        FUNCTYPE = ctypes.CFUNCTYPE(None, ctypes.c_char_p, ctypes.c_int, ctypes.c_int)

        self.lib.skip_final_au(
            ctypes.c_char_p(sequence.lower().encode()),
            ctypes.c_char_p(dot_bracket.encode()),
            ctypes.c_int(left_loop_stems),
            ctypes.c_int(right_loop_stems),
            FUNCTYPE(add_result),
        )

        return results
