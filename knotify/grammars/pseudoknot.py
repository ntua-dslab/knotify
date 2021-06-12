import itertools

import jinja2

TEMPLATE = """
S :
{% for idx, (p1, p2) in combinations %}
  | '{{ p1[0] }}' L '{{ p2[0] }}' D '{{ p1[1] }}' L '{{ p2[1] }}' # R{{ idx }} (1 3 5)
{% endfor %}
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
    allowed_pairs = {
        True: ["au", "gc", "ua", "cg", "gu", "ug"],
        False: ["au", "gc", "ua", "cg"],
    }[allow_ug]

    return (
        jinja2.Environment()
        .from_string(TEMPLATE)
        .render(
            combinations=enumerate(itertools.product(allowed_pairs, allowed_pairs)),
            max_dd_size=max_dd_size,
        )
    )
