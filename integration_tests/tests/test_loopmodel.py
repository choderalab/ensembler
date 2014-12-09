import os
from mock import Mock
import ensembler
import ensembler.initproject
import ensembler.integration_test_utils
import ensembler.tests.utils
import yaml


def test_find_loopmodel_executable():
    ensembler.core.find_loopmodel_executable()


def test_loopmodel_KC1D_HUMAN_D0_4KB8_D():
    with ensembler.integration_test_utils.integration_test_context(project_data='templates_resolved'):
        template = Mock()
        template.templateid = 'KC1D_HUMAN_D0_4KB8_D'
        template.chainid = 'D'
        template.resolved_seq = 'LRVGNRYRLGRKIGDIYLGTDIAAGEEVAIKLECPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDAQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNML'
        template.full_seq = 'LRVGNRYRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLKFGASRAADDAERERRDREERLRH'

        missing_residues = ensembler.initproject.pdbfix_template(template)
        ensembler.initproject.loopmodel_template(template, missing_residues)

        assert os.path.exists(os.path.join('templates', 'structures-modeled-loops', 'KC1D_HUMAN_D0_4KB8_D.pdb'))


def test_loopmodel_KC1D_HUMAN_D0_4HNF_A():
    with ensembler.integration_test_utils.integration_test_context(project_data='templates_resolved'):
        template = Mock()
        template.templateid = 'KC1D_HUMAN_D0_4HNF_A'
        template.chainid = 'A'
        template.resolved_seq = 'LRVGNRYRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLK'
        template.full_seq = 'MELRVGNRYRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLK'

        missing_residues = ensembler.initproject.pdbfix_template(template)
        ensembler.initproject.loopmodel_template(template, missing_residues)

        assert os.path.exists(os.path.join('templates', 'structures-modeled-loops', 'KC1D_HUMAN_D0_4HNF_A.pdb'))


def test_loopmodel_KC1D_HUMAN_D0_3UZP_A():
    """
    No missing residues
    """
    with ensembler.integration_test_utils.integration_test_context(project_data='templates_resolved'):
        template = Mock()
        template.templateid = 'KC1D_HUMAN_D0_3UZP_A'
        template.chainid = 'A'
        template.resolved_seq = 'LRVGNRYRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLK'
        template.full_seq = 'MELRVGNRYRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFRNLFHRQGFSYDYVFDWNMLK'

        missing_residues = ensembler.initproject.pdbfix_template(template)
        ensembler.initproject.loopmodel_template(template, missing_residues)

        assert not os.path.exists(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, 'KC1D_HUMAN_D0_3UZP_A.pdb'))
        log = yaml.load(open(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, 'KC1D_HUMAN_D0_3UZP_A-loopmodel-log.yaml')))
        assert log['no_missing_residues'] == True


def test_loopmodel_ZAP70_HUMAN_D0_1U59_A():
    with ensembler.integration_test_utils.integration_test_context(project_data='templates_resolved'):
        template = Mock()
        template.templateid = 'ZAP70_HUMAN_D0_1U59_A'
        template.chainid = 'A'
        template.resolved_seq = 'KKLFLKRDNLLIADIELGCGNFGSVRQGVYRKQIDVAIKVLKQGTEKADTEEMMREAQIMHQLDNPYIVRLIGVCQAEALMLVMEMAGGGPLHKFLVGKREEIPVSNVAELLHQVSMGMKYLEEKNFVHRDLAARNVLLVNRHYAKISDFGLSKALGADDSYYTARSAGKWPLKWYAPECINFRKFSSRSDVWSYGVTMWEALSYGQKPYKKMKGPEVMAFIEQGKRMECPPECPPELYALMSDCWIYKWEDRPDFLTVEQRMRACYYSLASKVEG'
        template.full_seq = 'DKKLFLKRDNLLIADIELGCGNFGSVRQGVYRMRKKQIDVAIKVLKQGTEKADTEEMMREAQIMHQLDNPYIVRLIGVCQAEALMLVMEMAGGGPLHKFLVGKREEIPVSNVAELLHQVSMGMKYLEEKNFVHRDLAARNVLLVNRHYAKISDFGLSKALGADDSYYTARSAGKWPLKWYAPECINFRKFSSRSDVWSYGVTMWEALSYGQKPYKKMKGPEVMAFIEQGKRMECPPECPPELYALMSDCWIYKWEDRPDFLTVEQRMRACYYSLASKVEG'

        missing_residues = ensembler.initproject.pdbfix_template(template)
        ensembler.initproject.loopmodel_template(template, missing_residues)

        assert os.path.exists(os.path.join('templates', 'structures-modeled-loops', 'ZAP70_HUMAN_D0_1U59_A.pdb'))