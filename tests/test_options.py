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
import os
from typing import Any

import pytest

from knotify import knotify
from knotify.algorithm.base import BaseAlgorithm


@pytest.mark.parametrize(
    "env_opt,env_val,opt,val",
    [
        ("KNOTIFY_ALLOW_SKIP_FINAL_AU", "True", "allow_skip_final_au", True),
        ("KNOTIFY_ALLOW_SKIP_FINAL_AU", "False", "allow_skip_final_au", False),
        ("KNOTIFY_ALLOW_UG", "1", "allow_ug", True),
        ("KNOTIFY_ALLOW_UG", "0", "allow_ug", False),
        ("KNOTIFY_MAX_HAIRPIN_BULGE", "2", "max_hairpin_bulge", 2),
        ("KNOTIFY_MAX_HAIRPIN_BULGE", "4", "max_hairpin_bulge", 4),
    ],
)
def test_environment(env_opt: str, env_val: str, opt: str, val: Any):
    opts = knotify.new_options()
    os.environ[env_opt] = env_val
    opts(args=[])

    algorithm, config = knotify.from_options(opts)

    assert isinstance(algorithm, BaseAlgorithm)
    assert config[opt] == val
