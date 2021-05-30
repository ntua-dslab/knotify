import jinja2

TEMPLATE = """
S : 'a' L 'a' D 'u' L 'u' # R1 (1 3 5)
  | 'u' L 'a' D 'a' L 'u' # R2 (1 3 5)
  | 'c' L 'a' D 'g' L 'u' # R3 (1 3 5)
  | 'g' L 'a' D 'c' L 'u' # R4 (1 3 5)
  | 'a' L 'u' D 'u' L 'a' # R7 (1 3 5)
  | 'u' L 'u' D 'a' L 'a' # R8 (1 3 5)
  | 'c' L 'u' D 'g' L 'a' # R9 (1 3 5)
  | 'g' L 'u' D 'c' L 'a' # R10 (1 3 5)
  | 'g' L 'u' D 'u' L 'a' # R11 (1 3 5)
  | 'a' L 'c' D 'u' L 'g' # R13 (1 3 5)
  | 'u' L 'c' D 'a' L 'g' # R14 (1 3 5)
  | 'c' L 'c' D 'g' L 'g' # R15 (1 3 5)
  | 'g' L 'c' D 'c' L 'g' # R16 (1 3 5)
  | 'a' L 'g' D 'u' L 'c' # R19 (1 3 5)
  | 'u' L 'g' D 'a' L 'c' # R20 (1 3 5)
  | 'c' L 'g' D 'g' L 'c' # R21 (1 3 5)
  | 'g' L 'g' D 'c' L 'c' # R22 (1 3 5)
{% if allow_ug %}
  | 'g' L 'a' D 'u' L 'u' # R5 (1 3 5)
  | 'u' L 'a' D 'g' L 'u' # R6 (1 3 5)
  | 'u' L 'u' D 'g' L 'a' # R12 (1 3 5)
  | 'g' L 'c' D 'u' L 'g' # R17 (1 3 5)
  | 'u' L 'c' D 'g' L 'g' # R18 (1 3 5)
  | 'u' L 'g' D 'g' L 'c' # R23 (1 3 5)
  | 'g' L 'g' D 'u' L 'c' # R24 (1 3 5)
  | 'a' L 'g' D 'u' L 'u' # R25 (1 3 5)
  | 'u' L 'g' D 'a' L 'u' # R26 (1 3 5)
  | 'c' L 'g' D 'g' L 'u' # R27 (1 3 5)
  | 'g' L 'g' D 'c' L 'u' # R28 (1 3 5)
  | 'u' L 'g' D 'g' L 'u' # R29 (1 3 5)
  | 'g' L 'g' D 'u' L 'u' # R30 (1 3 5)
  | 'a' L 'u' D 'u' L 'g' # R31 (1 3 5)
  | 'u' L 'u' D 'a' L 'g' # R32 (1 3 5)
  | 'c' L 'u' D 'g' L 'g' # R33 (1 3 5)
  | 'g' L 'u' D 'c' L 'g' # R34 (1 3 5)
{% endif %}
  ;

L : 'a' L # L1 (1)
  | 'u' L # L2 (1)
  | 'c' L # L3 (1)
  | 'g' L # L4 (1)
  | 'a'   # 0
  | 'u'   # 0
  | 'c'   # 0
  | 'g'   # 0
  ;

D : {% for x in range(max_dd_size) %}D{{ x }} {% endfor %} # M1 ({% for x in range(max_dd_size) %}{{ x }} {% endfor %})

E : 'a' # 0
  | 'u' # 0
  | 'c' # 0
  | 'g' # 0
  | ;

{% for x in range(max_dd_size) %}
D{{ x }} : E # N{{ x }} (0)
{% endfor %}
"""


def generate_grammar(allow_ug: bool, max_dd_size: int) -> str:
    """
    Generate grammar for pseudoknot detection, based on parameters.
    """
    return (
        jinja2.Environment()
        .from_string(TEMPLATE)
        .render(
            allow_ug=allow_ug,
            max_dd_size=max_dd_size,
        )
    )
