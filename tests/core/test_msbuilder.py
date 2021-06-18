from modelseedpy.core.msmodel import get_direction_from_constraints


def test_get_direction_from_constraints1():
    assert get_direction_from_constraints(-10, 10) == '='


def test_get_direction_from_constraints2():
    assert get_direction_from_constraints(0, 10) == '>'


def test_get_direction_from_constraints3():
    assert get_direction_from_constraints(5, 10) == '>'


def test_get_direction_from_constraints4():
    assert get_direction_from_constraints(-10, 0) == '<'


def test_get_direction_from_constraints5():
    assert get_direction_from_constraints(-10, -5) == '<'


def test_get_direction_from_constraints6():
    assert get_direction_from_constraints(0, 0) == '?'
