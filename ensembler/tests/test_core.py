import ensembler
import ensembler.param_parsers
import simtk.unit
from nose.plugins.attrib import attr


@attr('unit')
def test_project_metadata():
    ensembler.core.ProjectMetadata(project_stage='init')


@attr('unit')
def test_parse_api_params_string():
    params_dict = ensembler.param_parsers.parse_api_params_string('{"a": 3.2 / picoseconds, "b": "x", "c": 2.4}')
    assert params_dict == {'a': 3.2 / simtk.unit.picosecond, 'b': 'x', 'c': 2.4}


@attr('unit')
def test_eval_quantity_string():
    quantity = ensembler.param_parsers.eval_quantity_string('2 picoseconds')
    assert quantity == 2 * simtk.unit.picosecond