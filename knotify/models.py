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

# basic data structures

# flags tuple convenient tuple traversal --> defines the part
# of the tuple we need
OUTER = 1
INNER = 0


class Pknot(object):
    def __init__(
        self,
        len,
        left_core_indices=None,
        right_core_indices=None,
        left_loop_stems=None,
        right_loop_stems=None,
    ):
        """Define pseudoknot parts of a string
        Given a pseudoknot string we get the indexes that correspond to the
        structure bellow:
            "...LL.LLL(....RR.RR{..)LLL.LL....}RR.RR....."

        :param tuple left_core_indices: indexes of ()
        :param tuple right_core_indices: indexes of {}
        :param tuple left_loop_stems: indices of Ls ([inner], [outer])
        :param tuple right_loop_stems: indices of Rs ([inner], [outer])

        :param int len: the length of the pseudoknot
        """
        self._repr = {}
        if left_core_indices is None:
            self._repr["left_core_indices"] = []
        else:
            self._repr["left_core_indices"] = self._order_transparently(
                *left_core_indices, ascending=False
            )
        if right_core_indices is None:
            self._repr["right_core_indices"] = []
        else:
            self._repr["right_core_indices"] = self._order_transparently(
                *right_core_indices, ascending=True
            )
        if left_loop_stems is None:
            self._repr["left_loop_stems"] = []
        else:
            self._repr["left_loop_stems"] = self._order_transparently(
                *left_loop_stems, ascending=False
            )

        if right_loop_stems is None:
            self._repr["right_loop_stems"] = []
        else:
            self._repr["right_loop_stems"] = self._order_transparently(
                *right_loop_stems, ascending=True
            )

        self.len = len

    def serialize(self):
        """Returns the serialized vestion of the pseudoknot"""
        return self._repr

    @staticmethod
    def _order_transparently(el1, el2, ascending=None):
        # ascending should be boolean, None fallbacks to false
        # this is to make tuple assignments transparent
        s1 = el1
        s2 = el2

        # in case of None or empty Lists
        if not el1 and not el2:
            return (el1, el2)

        if type(el1) == list:
            s1 = el1[0]
            s2 = el2[0]

        if ascending:
            return (el1, el2) if s1 < s2 else (el2, el1)
        return (el1, el2) if s1 > s2 else (el2, el1)

    # core part
    @property
    def left_core_indices(self):
        """The left loop indices"""
        return self._repr.get("left_core_indices", None)

    @left_core_indices.setter
    def left_core_indices(self, bounds):
        inner, outer = bounds
        self._repr["left_core_indices"] = self._order_transparently(
            inner, outer, ascending=False
        )

    @property
    def right_core_indices(self):
        """The right loop indices"""
        return self._repr.get("right_core_indices", None)

    @right_core_indices.setter
    def right_core_indices(self, bounds):
        inner, outer = bounds
        self._repr["right_core_indices"] = self._order_transparently(
            inner, outer, ascending=True
        )

    # core stems
    @property
    def left_loop_stems(self):
        """
        :rtype: dict
        """
        return self._repr["left_loop_stems"]

    @left_loop_stems.setter
    def left_loop_stems(self, bounds):
        """
        :param list inner: the inner indices
        :param list outer: the outer indices
        """
        inner, outer = bounds
        self._repr["left_loop_stems"] = self._order_transparently(
            inner, outer, ascending=False
        )

    @property
    def right_loop_stems(self):
        """
        :rtype: dict
        """
        return self._repr["right_loop_stems"]

    @right_loop_stems.setter
    def right_loop_stems(self, bounds):
        """
        :param list inner: the inner indices
        :param list outer: the outer indices
        """
        inner, outer = bounds
        self._repr["right_loop_stems"] = self._order_transparently(
            inner, outer, ascending=True
        )

    # left over parts
    def get_left_outer_part(self):
        return list(
            set(range(0, self.left_core_indices[OUTER]))
            ^ set(self.left_loop_stems[OUTER])
        )

    def get_left_inner_part(self):
        return list(
            set(
                range(self.left_core_indices[OUTER] + 1, self.right_core_indices[INNER])
            )
            ^ set(self.right_loop_stems[INNER])
        )

    def get_right_outer_part(self):
        return list(
            set(range(self.right_core_indices[OUTER] + 1, self.len))
            ^ set(self.right_loop_stems[OUTER])
        )

    def get_right_inner_part(self):
        return list(
            set(
                range(self.left_core_indices[INNER] + 1, self.right_core_indices[OUTER])
            )
            ^ set(self.left_loop_stems[INNER])
        )

    def get_dd_seg(self):
        # FIXME: this does not work
        return range(self.right_core_indices[INNER] + 1, self.left_core_indices[INNER])

    # returns the indices of the pknot parts to check for basic stem alignments

    # FIXME: move that to a distincut utils functions
    def get_left_outer_potential(self):
        return (0, self.left_core_indices[OUTER])

    def get_right_inner_potential(self):
        return (self.left_core_indices[OUTER] + 1, self.right_core_indices[INNER])

    def get_right_outer_potential(self):
        return (self.right_core_indices[OUTER] + 1, self.len)

    def get_left_inner_potential(self):
        return (self.left_core_indices[INNER] + 1, self.right_core_indices[OUTER])

    # utility methods do follow:
    def get_left_loop_size(self):
        """Return the size of the left loop

        :rtype: int
        """
        return self.right_core_indices[INNER] - self.left_core_indices[OUTER] - 1

    def get_right_loop_size(self):
        """Return the size of the right loop

        :rtype: int
        """
        return self.right_core_indices[OUTER] - self.left_core_indices[INNER] - 1

    def to_dict(self):
        """
        Return pseudoknot information in dictionary.
        """

        return {
            "left_core_inner": self.left_core_indices[INNER],
            "left_core_outer": self.left_core_indices[OUTER],
            "right_core_inner": self.right_core_indices[INNER],
            "right_core_outer": self.right_core_indices[OUTER],
            "left_loop_size": self.get_left_loop_size(),
            "right_loop_size": self.get_right_loop_size(),
            "left_loop_stems": len(self.left_loop_stems[0]),
            "right_loop_stems": len(self.right_loop_stems[0]),
            "dot_bracket": self.visualize_dot_bracket(),
            "dd": str(len(self.get_dd_seg())),
        }

    def visualize_dot_bracket(self):
        """Serializes a pknot object to a dot_bracket representation

        :param Pknot pknot: the pknot to be serialized

        :returns: a dot bracket representation
        :rtype: str
        """
        dot_bracket = ["." for _ in range(0, self.len)]
        for i in self.left_loop_stems[INNER]:
            dot_bracket[i] = ")"
        for i in self.right_loop_stems[OUTER]:
            dot_bracket[i] = "]"
        for i in self.left_loop_stems[OUTER]:
            dot_bracket[i] = "("
        for i in self.right_loop_stems[INNER]:
            dot_bracket[i] = "["
        dot_bracket[self.left_core_indices[INNER]] = ")"
        dot_bracket[self.left_core_indices[OUTER]] = "("
        dot_bracket[self.right_core_indices[INNER]] = "["
        dot_bracket[self.right_core_indices[OUTER]] = "]"

        return "".join(dot_bracket)
