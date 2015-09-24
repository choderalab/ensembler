import os
import shutil
from simtk import unit
from ensembler.core import ProjectMetadata, manual_overrides_filename, ManualOverrides
from ensembler.core import get_valid_model_ids, default_project_dirnames
from ensembler.core import model_filenames_by_ensembler_stage
from ensembler.param_parsers import parse_api_params_string, eval_quantity_string
from ensembler.utils import enter_temp_dir, get_installed_resource_filename
from ensembler.tests.integrationtest_utils import integrationtest_context
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


@attr('unit')
def test_get_valid_model_filepaths():
    targetid = 'EGFR_HUMAN_D0'
    templateids = ['KC1D_HUMAN_D0_4KB8_D', 'KC1D_HUMAN_D0_4HNF_A']
    with integrationtest_context('refined_implicit'):
        valid_model_filenames = get_valid_model_ids('refine_implicit_md', targetid)
        assert all([fpath in templateids for fpath in valid_model_filenames])
