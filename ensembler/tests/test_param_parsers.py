from ensembler.param_parsers import eval_quantity_string, parse_api_params_string
from simtk import unit
from nose.plugins.attrib import attr


@attr('unit')
def test_eval_quantity_string():
    assert eval_quantity_string('2 picoseconds') == 2 * unit.picoseconds
    assert eval_quantity_string('2picoseconds') == 2 * unit.picoseconds
    assert eval_quantity_string('2 / picoseconds') == 2 / unit.picoseconds
    assert eval_quantity_string('2 nanosecond') == 2 * unit.nanoseconds


@attr('unit')
def test_parse_api_params_string():
    parsed = parse_api_params_string('{"a": 3 / picoseconds, "b": "x", "c": 2.4}')
    assert parsed == {
        'a': 3 / unit.picoseconds,
        'b': 'x',
        'c': 2.4
    }