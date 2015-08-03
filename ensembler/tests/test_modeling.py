import os
import shutil
import datetime
from mock import Mock
from nose.plugins.attrib import attr
import ensembler
import ensembler.tests
import ensembler.modeling
from ensembler.tests.utils import get_installed_resource_filename
from ensembler.tests.integrationtest_utils import integrationtest_context
import ensembler.cli_commands
from ensembler.utils import enter_temp_dir


@attr('modeller')
def test_import_modeller():
    import modeller
    import modeller.automodel


@attr('modeller')
def test_build_model():
    template_filepath = get_installed_resource_filename(os.path.join('resources',  'mock_template.pdb'))
    aln_filepath = get_installed_resource_filename(os.path.join('resources', 'mock_template-alignment.pir'))

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

        ensembler.modeling.build_model(target, template, target_setup_data=target_setup_data)

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


@attr('unit')
def test_align_command():
    ref_resources_dirpath = get_installed_resource_filename('example_project')
    with integrationtest_context(set_up_project_stage='templates_modeled_loops'):
        targets = ['KC1D_HUMAN_D0', 'EGFR_HUMAN_D0']
        templates = ['KC1D_HUMAN_D0_4KB8_D', 'KC1D_HUMAN_D0_4HNF_A']
        args = {
            '--targets': ','.join(targets),
            '--targetsfile': False,
            '--templates': ','.join(templates),
            '--templatesfile': False,
            '--verbose': False,
        }

        ensembler.cli_commands.align.dispatch(args)
        for target in targets:
            naln_files = 0
            for dir, subdirs, files in os.walk(os.path.join(ensembler.core.default_project_dirnames.models, target)):
                for file in files:
                    if file == 'alignment.pir':
                        naln_files += 1
            assert naln_files == len(templates)

        for target in targets:
            seqid_filepath = os.path.join(ensembler.core.default_project_dirnames.models, target, 'sequence-identities.txt')
            ref_seqid_filepath = os.path.join(ref_resources_dirpath, seqid_filepath)
            with open(seqid_filepath) as seqid_file:
                seqid_file_text = seqid_file.read()
            with open(ref_seqid_filepath) as ref_seqid_file:
                ref_seqid_file_text = ref_seqid_file.read()
            print(seqid_file_text)
            print(ref_seqid_file_text)
            assert seqid_file_text == ref_seqid_file_text


@attr('slow')
def test_build_models_command():
    with integrationtest_context(set_up_project_stage='aligned'):
        args = {
            '--targets': 'EGFR_HUMAN_D0',
            '--targetsfile': None,
            '--template_seqid_cutoff': None,
            '--templates': ','.join(['KC1D_HUMAN_D0_4KB8_D', 'KC1D_HUMAN_D0_4HNF_A']),
            '--templatesfile': None,
            '--write_modeller_restraints_file': None,
            '--verbose': False,
            '--help': False,
        }
        ensembler.cli_commands.build_models.dispatch(args)
        assert os.path.exists(os.path.join(ensembler.core.default_project_dirnames.models, 'EGFR_HUMAN_D0', 'KC1D_HUMAN_D0_4KB8_D', 'model.pdb'))
        assert os.path.exists(os.path.join(ensembler.core.default_project_dirnames.models, 'EGFR_HUMAN_D0', 'KC1D_HUMAN_D0_4HNF_A', 'model.pdb'))
        assert not os.path.exists(os.path.join(ensembler.core.default_project_dirnames.models, 'EGFR_HUMAN_D0', 'KC1D_HUMAN_D0_4KB8_D', 'restraints.rsr.gz'))
        assert not os.path.exists(os.path.join(ensembler.core.default_project_dirnames.models, 'EGFR_HUMAN_D0', 'KC1D_HUMAN_D0_4HNF_A', 'restraints.rsr.gz'))


@attr('unit')
def test_cluster_models():
    with integrationtest_context(set_up_project_stage='modeled'):
        ensembler.modeling.cluster_models()

@attr('unit')
def test_cluster_models_command():
    with integrationtest_context(set_up_project_stage='modeled'):
        args = {
            '--targetsfile': False,
            '--targets': False,
            '--cutoff': False,
            '--verbose': False,
            '--help': False,
        }
        ensembler.cli_commands.cluster.dispatch(args)