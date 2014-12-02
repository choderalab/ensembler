import msmseeder
import msmseeder.initproject
import gzip
import os
from mock import Mock

def test_pdbfix_KC1D_HUMAN_D0_4KB8_D():
    template_pdb_gz_filepath = os.path.abspath(os.path.join('tests', 'resources', 'KC1D_HUMAN_D0_4KB8_D.pdb.gz'))
    template_pdb_filepath = os.path.join(msmseeder.core.default_project_dirnames.templates_structures_resolved, 'KC1D_HUMAN_D0_4KB8_D.pdb')
    with msmseeder.utils.enter_temp_dir():
        msmseeder.utils.create_dir(msmseeder.core.default_project_dirnames.templates_structures_resolved)
        msmseeder.utils.create_dir(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops)
        with gzip.open(template_pdb_gz_filepath) as template_pdb_gz_file:
            with open(template_pdb_filepath, 'w') as template_pdb_file:
                template_pdb_file.write(template_pdb_gz_file.read())
        template = Mock()
        template.templateid = 'KC1D_HUMAN_D0_4KB8_D'
        template.chainid = 'D'
        template.resolved_seq = 'LRVGNRYRLGRKIGDIYLGTDIAAGEEVAIKLECPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDAQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNML'
        template.full_seq = 'LRVGNRYRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLKFGASRAADDAERERRDREERLRH'

        missing_residues = msmseeder.initproject.pdbfix_template(template)

        assert (0, 278) not in missing_residues
        assert missing_residues == {
            (0, 14): ['SER', 'GLY', 'SER', 'PHE', 'GLY'],
            (0, 34): ['VAL', 'LYS', 'THR', 'LYS', 'HIS'],
            (0, 147): ['ARG', 'THR', 'HIS'],
        }


def test_pdbfix_ABL1_HUMAN_D0_2E2B_B():
    template_pdb_gz_filepath = os.path.abspath(os.path.join('tests', 'resources', 'ABL1_HUMAN_D0_2E2B_B.pdb.gz'))
    template_pdb_filepath = os.path.join(msmseeder.core.default_project_dirnames.templates_structures_resolved, 'ABL1_HUMAN_D0_2E2B_B.pdb')
    with msmseeder.utils.enter_temp_dir():
        msmseeder.utils.create_dir(msmseeder.core.default_project_dirnames.templates_structures_resolved)
        msmseeder.utils.create_dir(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops)
        with gzip.open(template_pdb_gz_filepath) as template_pdb_gz_file:
            with open(template_pdb_filepath, 'w') as template_pdb_file:
                template_pdb_file.write(template_pdb_gz_file.read())
        template = Mock()
        template.templateid = 'ABL1_HUMAN_D0_2E2B_B'
        template.chainid = 'B'
        template.resolved_seq = 'ITMKHKLGGGQYGEVYEGVWKKYSLTVAVKTLEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIITEFMTYGNLLDYLRECNRQEVNAVVLLYMATQISSAMEYLEKKNFIHRDLAARNCLVGENHLVKVADFGLSTYTAHAGAKFPIKWTAPESLAYNKFSIKSDVWAFGVLLWEIATYGMSPYPGIDLSQVYELLEKDYRMERPEGCPEKVYELMRACWQWNPSDRPSFAEIHQAFETMF'
        template.full_seq = 'ITMKHKLGGGQYGEVYEGVWKKYSLTVAVKTLKEDTMEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIITEFMTYGNLLDYLRECNRQEVNAVVLLYMATQISSAMEYLEKKNFIHRDLAARNCLVGENHLVKVADFGLSRLMTGDTYTAHAGAKFPIKWTAPESLAYNKFSIKSDVWAFGVLLWEIATYGMSPYPGIDLSQVYELLEKDYRMERPEGCPEKVYELMRACWQWNPSDRPSFAEIHQAFETMFQESSISDEVEKELGKQ'

        missing_residues = msmseeder.initproject.pdbfix_template(template)

        assert (0, 271) not in missing_residues
        assert missing_residues == {
            (0, 32): ['LYS', 'GLU', 'ASP', 'THR', 'MET'],
            (0, 139): ['ARG', 'LEU', 'MET', 'THR', 'GLY', 'ASP'],
        }


def test_pdbfix_templates():
    template1_pdb_gz_filepath = os.path.abspath(os.path.join('tests', 'resources', 'KC1D_HUMAN_D0_4KB8_D.pdb.gz'))
    template1_pdb_filepath = os.path.join(msmseeder.core.default_project_dirnames.templates_structures_resolved, 'KC1D_HUMAN_D0_4KB8_D.pdb')
    template2_pdb_gz_filepath = os.path.abspath(os.path.join('tests', 'resources', 'KC1D_HUMAN_D0_3UYS_D.pdb.gz'))
    template2_pdb_filepath = os.path.join(msmseeder.core.default_project_dirnames.templates_structures_resolved, 'KC1D_HUMAN_D0_3UYS_D.pdb')
    with msmseeder.utils.enter_temp_dir():
        msmseeder.utils.create_dir(msmseeder.core.default_project_dirnames.templates_structures_resolved)
        msmseeder.utils.create_dir(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops)
        with gzip.open(template1_pdb_gz_filepath) as template1_pdb_gz_file:
            with open(template1_pdb_filepath, 'w') as template1_pdb_file:
                template1_pdb_file.write(template1_pdb_gz_file.read())
        with gzip.open(template2_pdb_gz_filepath) as template2_pdb_gz_file:
            with open(template2_pdb_filepath, 'w') as template2_pdb_file:
                template2_pdb_file.write(template2_pdb_gz_file.read())

        template1 = Mock()
        template1.templateid = 'KC1D_HUMAN_D0_4KB8_D'
        template1.chainid = 'D'
        template1.resolved_seq = 'LRVGNRYRLGRKIGDIYLGTDIAAGEEVAIKLECPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDAQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNML'
        template1.full_seq = 'LRVGNRYRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLKFGASRAADDAERERRDREERLRH'

        template2 = Mock()
        template2.templateid = 'KC1D_HUMAN_D0_3UYS_D'
        template2.chainid = 'D'
        template2.resolved_seq = 'RVGNRYRLGRKIGSDIYLGTDIAAGEEVAIKLECVKPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLK'
        template2.full_seq = 'MELRVGNRYRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLK'

        templates = [template1, template2]

        missing_residues_list = msmseeder.initproject.pdbfix_templates(templates)

        # assert (0, 278) not in missing_residues
        # assert missing_residues == {
        #     (0, 14): ['SER', 'GLY', 'SER', 'PHE', 'GLY'],
        #     (0, 34): ['VAL', 'LYS', 'THR', 'LYS', 'HIS'],
        #     (0, 147): ['ARG', 'THR', 'HIS'],
        # }


if __name__ == '__main__':
    test_pdbfix_templates()