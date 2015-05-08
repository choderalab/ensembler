import gzip
import os
import yaml
from mock import Mock
from nose.plugins.attrib import attr
import ensembler
import ensembler.initproject
from ensembler.tests.utils import get_installed_resource_filename
from ensembler.tests.integrationtest_utils import integrationtest_context
from ensembler.modeling import pdbfix_templates, pdbfix_template, loopmodel_template


@attr('unit')
def test_pdbfix_KC1D_HUMAN_D0_4KB8_D():
    template_pdb_gz_filepath = get_installed_resource_filename(os.path.join('resources', 'KC1D_HUMAN_D0_4KB8_D.pdb.gz'))
    template_pdb_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_resolved, 'KC1D_HUMAN_D0_4KB8_D.pdb')
    with ensembler.utils.enter_temp_dir():
        ensembler.utils.create_dir(ensembler.core.default_project_dirnames.templates_structures_resolved)
        ensembler.utils.create_dir(ensembler.core.default_project_dirnames.templates_structures_modeled_loops)
        with gzip.open(template_pdb_gz_filepath) as template_pdb_gz_file:
            with open(template_pdb_filepath, 'w') as template_pdb_file:
                contents = template_pdb_gz_file.read()
                if type(contents) == bytes:
                    contents = contents.decode('utf-8')
                template_pdb_file.write(contents)
        template = Mock()
        template.id = 'KC1D_HUMAN_D0_4KB8_D'
        template.seq = 'LRVGNRYRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLKFGASRAADDAERERRDREERLRH'

        missing_residues = pdbfix_template(template)

        assert (0, 278) not in missing_residues
        assert missing_residues == {
            (0, 14): ['SER', 'GLY', 'SER', 'PHE', 'GLY'],
            (0, 34): ['VAL', 'LYS', 'THR', 'LYS', 'HIS'],
            (0, 147): ['ARG', 'THR', 'HIS'],
        }


@attr('unit')
def test_pdbfix_ABL1_HUMAN_D0_2E2B_B():
    template_pdb_gz_filepath = get_installed_resource_filename(os.path.join('resources', 'ABL1_HUMAN_D0_2E2B_B.pdb.gz'))
    template_pdb_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_resolved, 'ABL1_HUMAN_D0_2E2B_B.pdb')
    with ensembler.utils.enter_temp_dir():
        ensembler.utils.create_dir(ensembler.core.default_project_dirnames.templates_structures_resolved)
        ensembler.utils.create_dir(ensembler.core.default_project_dirnames.templates_structures_modeled_loops)
        with gzip.open(template_pdb_gz_filepath) as template_pdb_gz_file:
            with open(template_pdb_filepath, 'w') as template_pdb_file:
                contents = template_pdb_gz_file.read()
                if type(contents) == bytes:
                    contents = contents.decode('utf-8')
                template_pdb_file.write(contents)
        template = Mock()
        template.id = 'ABL1_HUMAN_D0_2E2B_B'
        template.seq = 'ITMKHKLGGGQYGEVYEGVWKKYSLTVAVKTLKEDTMEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIITEFMTYGNLLDYLRECNRQEVNAVVLLYMATQISSAMEYLEKKNFIHRDLAARNCLVGENHLVKVADFGLSRLMTGDTYTAHAGAKFPIKWTAPESLAYNKFSIKSDVWAFGVLLWEIATYGMSPYPGIDLSQVYELLEKDYRMERPEGCPEKVYELMRACWQWNPSDRPSFAEIHQAFETMFQESSISDEVEKELGKQ'

        missing_residues = pdbfix_template(template)

        assert (0, 271) not in missing_residues
        assert missing_residues == {
            (0, 32): ['LYS', 'GLU', 'ASP', 'THR', 'MET'],
            (0, 139): ['ARG', 'LEU', 'MET', 'THR', 'GLY', 'ASP'],
        }


@attr('unit')
def test_pdbfix_templates():
    template1_pdb_gz_filepath = get_installed_resource_filename(os.path.join('resources', 'KC1D_HUMAN_D0_4KB8_D.pdb.gz'))
    template1_pdb_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_resolved, 'KC1D_HUMAN_D0_4KB8_D.pdb')
    template2_pdb_gz_filepath = get_installed_resource_filename(os.path.join('resources', 'KC1D_HUMAN_D0_3UYS_D.pdb.gz'))
    template2_pdb_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_resolved, 'KC1D_HUMAN_D0_3UYS_D.pdb')
    with ensembler.utils.enter_temp_dir():
        ensembler.utils.create_dir(ensembler.core.default_project_dirnames.templates_structures_resolved)
        ensembler.utils.create_dir(ensembler.core.default_project_dirnames.templates_structures_modeled_loops)
        with gzip.open(template1_pdb_gz_filepath) as template1_pdb_gz_file:
            with open(template1_pdb_filepath, 'w') as template1_pdb_file:
                contents = template1_pdb_gz_file.read()
                if type(contents) == bytes:
                    contents = contents.decode('utf-8')
                template1_pdb_file.write(contents)
        with gzip.open(template2_pdb_gz_filepath) as template2_pdb_gz_file:
            with open(template2_pdb_filepath, 'w') as template2_pdb_file:
                contents = template2_pdb_gz_file.read()
                if type(contents) == bytes:
                    contents = contents.decode('utf-8')
                template2_pdb_file.write(contents)

        template1 = Mock()
        template1.id = 'KC1D_HUMAN_D0_4KB8_D'
        template1.seq = 'LRVGNRYRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLKFGASRAADDAERERRDREERLRH'

        template2 = Mock()
        template2.id = 'KC1D_HUMAN_D0_3UYS_D'
        template2.seq = 'MELRVGNRYRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLK'

        templates = [template1, template2]

        missing_residues_list = pdbfix_templates(templates)

        # assert (0, 278) not in missing_residues
        # assert missing_residues == {
        #     (0, 14): ['SER', 'GLY', 'SER', 'PHE', 'GLY'],
        #     (0, 34): ['VAL', 'LYS', 'THR', 'LYS', 'HIS'],
        #     (0, 147): ['ARG', 'THR', 'HIS'],
        # }


@attr('non_conda_dependencies')
def test_find_loopmodel_executable():
    ensembler.core.find_loopmodel_executable()


@attr('slow')
@attr('non_conda_dependencies')
def test_loopmodel_KC1D_HUMAN_D0_4KB8_D():
    with integrationtest_context(set_up_project_stage='templates_resolved'):
        template = Mock()
        template.id = 'KC1D_HUMAN_D0_4KB8_D'
        template.seq = 'YRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLF'

        missing_residues = pdbfix_template(template)
        loopmodel_template(template, missing_residues)

        assert os.path.exists(os.path.join('templates', 'structures-modeled-loops', 'KC1D_HUMAN_D0_4KB8_D.pdb'))


@attr('slow')
@attr('non_conda_dependencies')
def test_loopmodel_KC1D_HUMAN_D0_4HNF_A():
    with integrationtest_context(set_up_project_stage='templates_resolved'):
        template = Mock()
        template.id = 'KC1D_HUMAN_D0_4HNF_A'
        template.seq = 'YRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLF'

        missing_residues = pdbfix_template(template)
        loopmodel_template(template, missing_residues)

        assert os.path.exists(os.path.join('templates', 'structures-modeled-loops', 'KC1D_HUMAN_D0_4HNF_A.pdb'))


@attr('slow')
@attr('non_conda_dependencies')
def test_loopmodel_KC1D_HUMAN_D0_3UZP_A():
    """
    No missing residues
    """
    with integrationtest_context(set_up_project_stage='templates_resolved'):
        template = Mock()
        template.id = 'KC1D_HUMAN_D0_3UZP_A'
        template.seq = 'YRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLF'

        missing_residues = pdbfix_template(template)
        loopmodel_template(template, missing_residues)

        assert not os.path.exists(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, 'KC1D_HUMAN_D0_3UZP_A.pdb'))
        log = yaml.load(open(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, 'KC1D_HUMAN_D0_3UZP_A-loopmodel-log.yaml')))
        assert log['no_missing_residues'] == True