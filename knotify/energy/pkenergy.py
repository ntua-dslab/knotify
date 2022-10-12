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

from knotify.energy.base import BaseEnergy


class PKEnergy(BaseEnergy):
    """
    Load MFE calculator from a dynamic library. The library should export:

    ```
    // will be called once during instantiation of the class. use to load any
    // configuation files (e.g. parameters).
    void initialize(char *config_dir, char *model);

    // will be called for each sequence.
    float get_energy(char *sequence, char *structure);
    ```
    """

    def __init__(self, library: str, config_dir: str, model: str):
        self._lib = ctypes.CDLL(library)
        self._lib.get_energy.restype = ctypes.c_double
        self._lib.initialize(
            ctypes.c_char_p(config_dir.encode()),
            ctypes.c_char_p(model.encode()),
        )

    def eval(self, sequence: str, dot_bracket: str) -> float:
        return self._lib.get_energy(
            ctypes.c_char_p(sequence.encode()),
            ctypes.c_char_p(dot_bracket.encode()),
        )
