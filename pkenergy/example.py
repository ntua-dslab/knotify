#!/usr/bin/env python3
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

c = ctypes.CDLL("./hotknots/LE/libpkenergy.so")

c.get_energy.restype = ctypes.c_double


def cstr(s: str):
    return ctypes.c_char_p(s.encode())


c.initialize(cstr("hotknots/params", cstr("dp")))

for x in range(1000):
    print(
        c.get_energy(
            cstr(
                "AAUGCAACUUUUAAAUAGUUUAUCUGUUAAGAUAAACCACCUAGGUUGCAUAUAUAAAAAAUAAAAGGUGCC"
            ),
            cstr(
                ".(((((((((...........................[[[[[)))))))))..............]]]]].."
            ),
        )
    )
