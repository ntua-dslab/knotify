from typing import Dict, Tuple
from dataclasses import dataclass

import pytest

from pseudoknot_detector.scoring import (
    Dot2Pair,
    find_matches,
    get_confusion_matrix,
    get_core_stem_indices,
    get_correct_core_stems,
)


@dataclass
class Case:
    dot_bracket: str
    result: Dict[str, int] = None
    core_stems: Tuple[int, int, int, int] = None
    raises: Exception = None


@pytest.mark.parametrize(
    "name,case",
    {
        "only_left_stems": Case(
            dot_bracket=".(.).",
            result={
                "right_core_stems": [],
                "left_core_stems": [],
                "right_stems": [],
                "left_stems": [(1, 3)],
            },
        ),
        "only_core_stems": Case(
            dot_bracket=".([)].",
            result={
                "right_core_stems": [(2, 4)],
                "left_core_stems": [(1, 3)],
                "right_stems": [],
                "left_stems": [],
            },
            core_stems=(1, 2, 3, 4),
        ),
        "core_and_left_stems": Case(
            dot_bracket="(([))]",
            result={
                "right_core_stems": [(2, 5)],
                "left_core_stems": [(1, 3)],
                "right_stems": [],
                "left_stems": [(0, 4)],
            },
            core_stems=(1, 2, 3, 5),
        ),
        "complete": Case(
            dot_bracket="(([[))]]",
            result={
                "right_core_stems": [(3, 6)],
                "left_core_stems": [(1, 4)],
                "right_stems": [(2, 7)],
                "left_stems": [(0, 5)],
            },
            core_stems=(1, 3, 4, 6),
        ),
        "alternate_symbols": Case(
            dot_bracket="[[((]]))",
            result={
                "right_core_stems": [(3, 6)],
                "left_core_stems": [(1, 4)],
                "right_stems": [(2, 7)],
                "left_stems": [(0, 5)],
            },
            core_stems=(1, 3, 4, 6),
        ),
    }.items(),
)
def test_parse_dot(name: str, case: Case):
    parser = Dot2Pair(case.dot_bracket)
    parser.parse_dot()

    assert parser.result == case.result

    if case.core_stems is not None:
        assert get_core_stem_indices(case.dot_bracket) == case.core_stems


@pytest.mark.parametrize(
    "name,case",
    {
        "no_stems": Case(
            dot_bracket=".......",
            raises=ValueError,
        ),
        "right_open_no_close": Case(
            dot_bracket="[[(((]]))",
            raises=ValueError,
        ),
        "right_close_no_open": Case(
            dot_bracket="[[((]])))",
            raises=IndexError,
        ),
        "left_open_no_close": Case(
            dot_bracket="[[[[((]]))",
            raises=ValueError,
        ),
        "left_close_no_open": Case(
            dot_bracket="[[((]]]]))",
            raises=IndexError,
        ),
        "left_close_no_open_with_neutral": Case(
            dot_bracket="[.........[((]]::::]]....))::::",
            raises=IndexError,
        ),
    }.items(),
)
def test_exceptions(name: str, case: Case):
    with pytest.raises(case.raises):
        parser = Dot2Pair(case.dot_bracket)
        parser.parse_dot()


@pytest.mark.parametrize(
    "bracket, result",
    [
        (
            "......",
            {},
        ),
        (
            "[[((]]))",
            {
                0: 5,
                1: 4,
                2: 7,
                3: 6,
                4: 1,
                5: 0,
                6: 3,
                7: 2,
            },
        ),
        (
            "[[((...]]))",
            {
                0: 8,
                1: 7,
                2: 10,
                3: 9,
                7: 1,
                8: 0,
                9: 3,
                10: 2,
            },
        ),
    ],
)
def test_pseudoknot(bracket, result):
    assert find_matches(bracket) == result


@pytest.mark.parametrize(
    "truth, pred, tp, tn, fp, fn, slack",
    [
        ("[(])", "[(])", 4, 0, 0, 0, 0),
        ("[(...])", "[(...])", 4, 3, 0, 0, 0),
        (".((((..[[[)))).....]]]..", ".((((..[[[)))).....]]]..", 14, 10, 0, 0, 0),
        ("..(((..[[[)))......]]]..", ".((((..[[[)))).....]]]..", 12, 10, 2, 0, 0),
        (".((((..[[[)))).....]]]..", "..(((..[[[)))......]]]..", 12, 10, 0, 2, 0),
        (
            ".............(((((..[[[.))))).........]]]..",
            "............((((((..[[[.))))))........]]]..",
            16,
            25,
            2,
            0,
            0,
        ),
        (
            ".............(((((..[[[.))))).........]]]..",
            ".............(((((..[[[.))))).........]]]..",
            16,
            27,
            0,
            0,
            0,
        ),
        # (
        #     # Disabled case: see TODO in get_confusion_matrix
        #     ".[(])",
        #     "[.(])",
        #     4,
        #     1,
        #     0,
        #     0,
        #     1,
        # ),
    ],
)
def test_confusion_matrix(
    truth: str, pred: str, tp: int, tn: int, fp: int, fn: int, slack: int
):
    assert get_confusion_matrix(truth, pred, slack=slack) == (tp, tn, fp, fn)


@pytest.mark.parametrize(
    "truth, pred, slack, result",
    [
        (
            "([)]",
            "([)]",
            0,
            2,
        ),
        (
            ".([)]",
            "(.[)]",
            0,
            1,
        ),
        (
            ".([)]",
            "(.[)]",
            1,
            2,
        ),
        (
            ".([).]",
            "(.[).]",
            1,
            2,
        ),
        (
            ".([).]",
            "(.[.)]",
            1,
            1,
        ),
        (
            ".([).]",
            "(.[.)]",
            0,
            1,
        ),
        (
            ".([).]..",
            "(.[.)..]",
            1,
            0,
        ),
        (
            "...([).]..",
            "..([..).].",
            1,
            0,
        ),
    ],
)
def test_correct_core_stems(truth: str, pred: str, slack: int, result: int):
    assert get_correct_core_stems(truth, pred, slack) == result
