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
class BaseParser:
    """
    Base RNA parser class. Defines the interface for using the RNA parser, and
    common configuration options for all parsers. Parser implementations are expected
    to read and respect these configurations.

    Implementations should
    """

    def __init__(
        self,
        max_dd_size: int = 2,
        allow_ug: bool = False,
        min_dd_size: int = 0,
        max_window_size: int = 100,
        min_window_size: int = 6,
        max_window_size_ratio: float = 0,
        min_window_size_ratio: float = 0,
    ):
        self.allow_ug = allow_ug
        self.max_dd_size = max_dd_size
        self.min_dd_size = min_dd_size
        self.max_window_size = max_window_size
        self.min_window_size = min_window_size
        self.max_window_size_ratio = max_window_size_ratio
        self.min_window_size_ratio = min_window_size_ratio

    def detect_pseudoknots(self, sequence: str) -> list:
        """
        Return a list of pseudoknots for the given RNA sequence (and any subsequences).
        The format should be:
        [
            "<left1>,<size1>,<left_loop_size1>,<dd_size1>",
            "<left2>,<size2>,<left_loop_size2>,<dd_size2>",
            ...
        ]
        """
        raise NotImplementedError
