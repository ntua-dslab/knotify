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

from knotify.parsers.base import BaseParser


class CTypesParser(BaseParser):
    """
    Base class for ctypes-based parsers. The library_path argument should be
    a shared library exposing the following functions with C linkage:

    ```
    // Will be called once to initialize the library. options is a string
    // with implementation-specific configuration.
    void initialize(
        char *options,
        int allow_ug,
        int min_dd_size,
        int max_dd_size,
        int min_window_size,
        int max_window_size
    );

    // Parse RNA sequence and return positions of possible core stems pseudoknots.
    // Returned as a char *buffer with the following format:
    //
    // ```
    // <left1>,<size1>,<leftloopsize1>,<ddsize1>
    // <left2>,<size2>,<leftloopsize2>,<ddsize2>
    // ...
    // ```
    char *detect_pseudoknots(char *sequence);
    ```
    """

    def get_options(self) -> str:
        raise NotImplementedError

    def __init__(self, library_path: str, *args, **kwargs):
        super(CTypesParser, self).__init__(*args, **kwargs)

        self.lib = ctypes.CDLL(library_path)
        self.lib.detect_pseudoknots.restype = ctypes.c_char_p
        self.lib.initialize(
            ctypes.c_char_p(self.get_options().encode()),
            ctypes.c_int(self.allow_ug),
            ctypes.c_int(self.min_dd_size),
            ctypes.c_int(self.max_dd_size),
            ctypes.c_int(self.min_window_size),
            ctypes.c_int(self.max_window_size),
            ctypes.c_float(self.min_window_size_ratio),
            ctypes.c_float(self.max_window_size_ratio),
        )

    def detect_pseudoknots(self, sequence: str) -> list:
        return (
            self.lib.detect_pseudoknots(ctypes.c_char_p(sequence.lower().encode()))
            .decode()
            .rstrip("\n")
            .split("\n")
        )
