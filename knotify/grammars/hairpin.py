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
