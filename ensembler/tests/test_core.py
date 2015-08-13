import os
import shutil
from ensembler.core import ProjectMetadata, manual_overrides_filename, ManualOverrides
from ensembler.param_parsers import parse_api_params_string, eval_quantity_string
from simtk import unit
from ensembler.utils import enter_temp_dir, get_installed_resource_filename
from nose.plugins.attrib import attr


@attr('unit')
def test_import_ensembler_version():
    import ensembler.version


@attr('unit')
def test_project_metadata():
    ProjectMetadata(project_stage='init')


@attr('unit')
def test_parse_api_params_string():
    params_dict = parse_api_params_string('{"a": 3.2 / picoseconds, "b": "x", "c": 2.4}')
    assert params_dict == {'a': 3.2 / unit.picosecond, 'b': 'x', 'c': 2.4}


@attr('unit')
def test_eval_quantity_string():
    quantity = eval_quantity_string('2 picoseconds')
    assert quantity == 2 * unit.picosecond


@attr('unit')
def test_manual_overrides_file():
    with enter_temp_dir():
        ref_manual_overrides_file = get_installed_resource_filename(
            os.path.join('tests', 'example_project', 'manual-overrides.yaml')
        )
        shutil.copy(ref_manual_overrides_file, manual_overrides_filename)
        manual_overrides = ManualOverrides()
        assert manual_overrides.target.domain_spans == {'ABL1_HUMAN_D0': '242-513'}
        assert manual_overrides.template.min_domain_len == 0
        assert manual_overrides.template.max_domain_len == 350
        assert manual_overrides.template.domain_spans == {'ABL1_HUMAN_D0': '242-513'}
        assert manual_overrides.template.skip_pdbs == [
            '4CYJ', '4P41', '4P2W', '4QTD', '4Q2A', '4CTB', '4QOX'
        ]
        assert manual_overrides.refinement.ph == 8.0
        assert manual_overrides.refinement.custom_residue_variants_by_targetid == {
            'EGFR_HUMAN_D0': {49: 'ASH'}
        }
