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

from knotify.pairalign.base import BasePairAlign


class CTypesPairAlign(BasePairAlign):
    """
    Consecutive RNA pairalign class. Match as many consecutive loop stems as possible.

    Implementations should be done in C code, exporting the following symbol:

    ```c
    void pairalign(char* sequence, int i, int j, int left_loop_size, int dd_size,
                   void (*cb)(char*, int, int));
    ```

    The implementation is done in C code in pairalign/cpairalign.c

    For usage, refer to the unit tests in test/test_pairalign.py
    """

    def __init__(self, library_path: str, *args, **kwargs):
        super(CTypesPairAlign, self).__init__(*args, **kwargs)

        self.lib = ctypes.CDLL(library_path)

    def pairalign(
        self, sequence: str, i: int, j: int, left_loop_size: int, dd_size: int
    ) -> str:
        results = []

        def add_result(dot_bracket, left_loop_stems, right_loop_stems):
            results.append((dot_bracket.decode(), left_loop_stems, right_loop_stems))

        FUNCTYPE = ctypes.CFUNCTYPE(None, ctypes.c_char_p, ctypes.c_int, ctypes.c_int)

        self.lib.pairalign(
            ctypes.c_char_p(sequence.lower().encode()),
            ctypes.c_int(i),
            ctypes.c_int(j),
            ctypes.c_int(left_loop_size),
            ctypes.c_int(dd_size),
            FUNCTYPE(add_result),
        )

        return results
