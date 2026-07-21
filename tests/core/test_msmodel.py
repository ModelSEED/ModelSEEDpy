# -*- coding: utf-8 -*-
from modelseedpy.core.msmodel import *


def test_get_direction_from_constraints1():
    res = get_direction_from_constraints(0, 1000)

    assert res == ">"


def test_get_direction_from_constraints2():
    res = get_direction_from_constraints(-1000, 0)

    assert res == "<"


def test_get_direction_from_constraints3():
    res = get_direction_from_constraints(-1000, 1000)

    assert res == "="


def test_get_set_set1():
    res = get_set_set("A")

    assert len(res) == 1
    assert {"A"} in res


def test_get_set_set2():
    res = get_set_set("A and B")

    assert len(res) == 1
    assert {"A", "B"} in res


def test_get_set_set3():
    res = get_set_set("A or B")

    assert len(res) == 2
    assert {"A"} in res
    assert {"B"} in res


def test_get_set_set4():
    res = get_set_set("A or B or C")

    assert len(res) == 3
    assert {"A"} in res
    assert {"B"} in res
    assert {"C"} in res


def test_get_set_set5():
    res = get_set_set("A or B and C")

    assert len(res) == 2
    assert {"A"} in res
    assert {"B", "C"} in res


def test_get_set_set6():
    res = get_set_set("A and B or C")

    assert len(res) == 2
    assert {"A", "B"} in res
    assert {"C"} in res


def test_get_set_set7():
    res = get_set_set("(A or B) and C")

    assert len(res) == 2
    assert {"A", "C"} in res
    assert {"B", "C"} in res


def test_get_set_set8():
    res = get_set_set("A and (B or C)")

    assert len(res) == 2
    assert {"A", "B"} in res
    assert {"A", "C"} in res
