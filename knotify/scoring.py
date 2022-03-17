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
from typing import Tuple, Dict
from enum import Enum


class State(Enum):
    NONE = 0
    LEFT_CORE = 1
    RIGHT_CORE = 2


class Stem(Enum):
    RIGHT = "right"
    LEFT = "left"


class Dot2Pair:
    """
    Parse dot bracket expressions and extract positions of stems
    """

    SYMBOL_MAP = {"(": ")", "[": "]", ")": "(", "]": "["}
    OPPOSITE_SYMBOL_MAP = {")": "]", "(": "[", "[": "(", "]": ")"}

    def __init__(self, dot_bracket: str):
        """A new parser set to parse the related dot bracket.

        :return:  a DotParser instance
        """
        self.dot_bracket = dot_bracket
        # these flags are going to be set
        # whenever we have a match a pseudoknot lookup ...
        self.state = State.NONE
        self.stack = {Stem.RIGHT: [], Stem.LEFT: []}

        try:
            first_symbol = next(x for x in self.dot_bracket if x in self.SYMBOL_MAP)
        except StopIteration:
            raise ValueError(f"{self.dot_bracket} has no stem symbols")

        self.left_open_symbol = first_symbol
        self.right_open_symbol = self.OPPOSITE_SYMBOL_MAP[first_symbol]
        self.left_close_symbol = self.SYMBOL_MAP[self.left_open_symbol]
        self.right_close_symbol = self.SYMBOL_MAP[self.right_open_symbol]

        # result stack
        self._result = {
            "right_core_stems": [],
            "left_core_stems": [],
            "right_stems": [],
            "left_stems": [],
        }

    def parse_dot(self):
        for position, char in enumerate(self.dot_bracket):
            self.__char_digest(position, char)

        if any(self.stack.values()):
            raise ValueError(f"imbalanced stems {self.stack}")

    def __char_digest(self, position: int, char: str):
        """
        :return: side-effects
        """
        if char == self.left_open_symbol:
            self.stack[Stem.LEFT].append(position)

        elif char == self.left_close_symbol:
            binding = self.stack[Stem.LEFT].pop()
            is_core = self.state == State.LEFT_CORE
            if is_core:
                self.state = State.RIGHT_CORE

            self.__add_pair(Stem.LEFT, (binding, position), is_core)

        elif char == self.right_open_symbol:
            self.stack[Stem.RIGHT].append(position)
            if self.state == State.NONE:
                self.state = State.LEFT_CORE

        elif char == self.right_close_symbol:
            binding = self.stack[Stem.RIGHT].pop()
            is_core = self.state == State.RIGHT_CORE
            if is_core:
                self.state = State.NONE

            self.__add_pair(Stem.RIGHT, (binding, position), is_core)

    def __add_pair(self, key: Stem, pair: Tuple[int, int], is_core: bool):
        key = f"{key.value}{'_core' if is_core else ''}_stems"
        self._result[key].append(pair)

    @property
    def result(self):
        return self._result


def find_matches(dot_bracket: str) -> Dict[int, int]:
    """
    parse dot bracket. return dictionary of matching stems. If no stems exist,
    an empty dictionary is returned.
    """
    try:
        dot2pair = Dot2Pair(dot_bracket)
        dot2pair.parse_dot()
    except ValueError:
        return {}

    matches = {}
    for stems in dot2pair.result.values():
        for start, end in stems:
            matches[start] = end
            matches[end] = start

    return matches


def get_confusion_matrix(truth: str, prediction: str, slack: int = 0):
    truth_matches = find_matches(truth)
    prediction_matches = find_matches(prediction)

    assert len(truth) == len(prediction)

    tp, tn, fp, fn = 0, 0, 0, 0

    for idx in range(0, len(truth)):
        tmatch = truth_matches.get(idx, None)
        pmatch = prediction_matches.get(idx, None)

        if tmatch is None:
            if pmatch is None:
                tn += 1
            else:
                # TODO: check disabled test, at index 0
                fp += 1
        else:
            if pmatch is None:
                fn += 1
                # TODO: check disabled test, at index 1
            else:
                if abs(pmatch - tmatch) <= slack:
                    tp += 1
                else:
                    fp += 1

    return tp, tn, fp, fn


def get_core_stem_indices(dot_bracket: str):
    """
    Returns core stem indices for the dot bracket. If pred does not
    have a pseudoknot, a ValueError or an IndexError is raised.
    """
    dot2pair = Dot2Pair(dot_bracket)
    dot2pair.parse_dot()

    result = dot2pair.result
    right_core = result["right_core_stems"][0]
    left_core = result["left_core_stems"][0]

    return left_core[0], right_core[0], left_core[1], right_core[1]


def get_correct_core_stems(truth: str, pred: str, slack: int = 0) -> int:
    """
    Returns number of correct core stems for the prediction. If there is
    pseudoknot in the predicted dot bracket, a ValueError is raised.
    """
    try:
        tstems = get_core_stem_indices(truth)
        pstems = get_core_stem_indices(pred)
    except IndexError:
        return 0

    left_ok = abs(tstems[0] - pstems[0]) + abs(tstems[2] - pstems[2]) <= slack
    right_ok = abs(tstems[1] - pstems[1]) + abs(tstems[3] - pstems[3]) <= slack

    return left_ok + right_ok
