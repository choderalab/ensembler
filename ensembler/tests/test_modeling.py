import os
import shutil
import datetime

from mock import Mock
from nose.plugins.attrib import attr

import ensembler
import ensembler.tests
import ensembler.modeling
import tests.integration_test_utils
import ensembler.cli_commands
from ensembler.utils import enter_temp_dir


@attr('unit')
def test_build_model():
    template_filepath = os.path.abspath(os.path.join('tests', 'resources', 'mock_template.pdb'))
    aln_filepath = os.path.abspath(os.path.join('tests', 'resources', 'mock_template-alignment.pir'))

    with enter_temp_dir():
        target = Mock()
        template = Mock()
        target_setup_data = Mock()
        target.id = 'mock_target'
        target.seq = 'YILGDTLGVGGKVKVGKH'
        template.id = 'mock_template'
        template.seq = 'YQNLSPVGSGGSVCAAFD'
        target_setup_data.target_starttime = datetime.datetime.utcnow()
        target_setup_data.models_target_dir = os.path.join(ensembler.core.default_project_dirnames.models, target.id)
        os.mkdir(ensembler.core.default_project_dirnames.models)
        model_dir = os.path.join(ensembler.core.default_project_dirnames.models, 'mock_target', 'mock_template')
        os.makedirs(model_dir)
        os.mkdir(ensembler.core.default_project_dirnames.templates)
        os.mkdir(ensembler.core.default_project_dirnames.templates_structures_resolved)

        shutil.copy(template_filepath, ensembler.core.default_project_dirnames.templates_structures_resolved)
        shutil.copy(aln_filepath, os.path.join(model_dir, 'alignment.pir'))

        ensembler.modeling.build_model(target, template, template, target_setup_data=target_setup_data)

        # Example model.pdb.gz contents (not testing this as it may be dependent
        # upon Modeller version as well as modelling stochasticity):
        #
        # ..EXPDTA    THEORETICAL MODEL, MODELLER 9.12 2014/08/26 13:15:44
        # REMARK   6 MODELLER OBJECTIVE FUNCTION:       326.6798
        # REMARK   6 MODELLER BEST TEMPLATE % SEQ ID:  27.778
        # ATOM      1  N   TYR     1      48.812  50.583  13.949  1.00110.28           N
        # ATOM      2  CA  TYR     1      49.070  50.334  15.387  1.00110.28           C

        model_filepath = os.path.join(target_setup_data.models_target_dir, template.id, 'model.pdb.gz')
        assert os.path.exists(model_filepath)
        assert os.path.getsize(model_filepath) > 0


@attr('integration')
def test_align_command():
    ref_resources_dirpath = os.path.abspath(os.path.join('tests', 'integration_test_resources'))
    with tests.integration_test_utils.integration_test_context(set_up_project_stage='templates_modeled_loops'):
        targets = ['KC1D_HUMAN_D0', 'EGFR_HUMAN_D0']
        args = {
            '--targets': targets,
            '--verbose': False,
        }
        ensembler.cli_commands.align.dispatch(args)
        for target in targets:
            naln_files = 0
            for dir, subdirs, files in os.walk(os.path.join(ensembler.core.default_project_dirnames.models, target)):
                for file in files:
                    if file == 'alignment.pir':
                        naln_files += 1
            assert naln_files == 67

        for target in targets:
            seqid_filepath = os.path.join(ensembler.core.default_project_dirnames.models, target, 'sequence-identities.txt')
            ref_seqid_filepath = os.path.join(ref_resources_dirpath, seqid_filepath)
            with open(seqid_filepath) as seqid_file:
                seqid_file_text = seqid_file.read()
            with open(ref_seqid_filepath) as ref_seqid_file:
                ref_seqid_file_text = ref_seqid_file.read()
            assert seqid_file_text == ref_seqid_file_text


@attr('integration')
def test_build_models_command():
    with tests.integration_test_utils.integration_test_context(set_up_project_stage='templates_resolved'):
        args = {
            '--targets': ['SRC_HUMAN_D0'],
            '--templates': ['KC1D_HUMAN_D0_3UZP_A', 'KC1D_HUMAN_D0_4KB8_D'],
            '--verbose': False,
            '--help': False,
        }
        ensembler.cli_commands.build_models.dispatch(args)
        assert os.path.exists(os.path.join(ensembler.core.default_project_dirnames.models, 'SRC_HUMAN_D0', 'KC1D_HUMAN_D0_3UZP_A', 'model.pdb'))
        assert os.path.exists(os.path.join(ensembler.core.default_project_dirnames.models, 'SRC_HUMAN_D0', 'KC1D_HUMAN_D0_4KB8_D', 'model.pdb'))