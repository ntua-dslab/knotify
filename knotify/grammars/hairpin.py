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
import jinja2

# Description for grammar:
#
# .(((....))).
# LPPPMMMMPPPR
#
# Initial something:
# - @start: index of first P
# - @stems: number of Ps (total of hairpin div 2)
# - @size: number of Ms
#
# For example above, SHOULD be: start=1, stems=3, size=4
#
# The reference implementation of the grammar tree parsing can be found
# in parsers/hairpin.c

TEMPLATE = """

S : L P R # S1 (0 1 2)

P : 'a' P 'u' # P1 (1)
  | 'u' P 'a' # P2 (1)
  | 'g' P 'c' # P3 (1)
  | 'c' P 'g' # P4 (1)
{% if allow_ug %}
  | 'g' P 'u' # P5 (1)
  | 'u' P 'g' # P6 (1)
{% endif %}
  | M # P5 (0)

L : K # L1 (0)

R : K # R1 (0)

M : K # M1 (0)

K : 'a' K # K1 (1)
  | 'u' K # K2 (1)
  | 'c' K # K3 (1)
  | 'g' K # K4 (1)
  | ;
"""


def generate_grammar(allow_ug: bool) -> str:
    """
    Generate grammar for hairpin detection, based on parameters.
    """
    return (
        jinja2.Environment()
        .from_string(TEMPLATE)
        .render(
            allow_ug=allow_ug,
        )
    )
