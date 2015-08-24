import os
from ensembler.tools.inspect import LoopmodelLogs
from ensembler.tools.mktraj import MkTrajImplicitStart
from ensembler.tools.quick_model import QuickModel
from ensembler.core import default_project_dirnames
from simtk import unit
from ensembler.utils import enter_temp_dir
from ensembler.tests.integrationtest_utils import integrationtest_context
from nose.plugins.attrib import attr


@attr('unit')
def test_loopmodel_logs():
    with integrationtest_context(set_up_project_stage='templates_modeled_loops'):
        loopmodel_logs = LoopmodelLogs()
        loopmodel_logs.add_missing_resis_data()


@attr('unit')
def test_mktraj_implicit_start():
    with integrationtest_context(set_up_project_stage='refined_explicit'):
        MkTrajImplicitStart(targetid='EGFR_HUMAN_D0', loglevel='debug')


@attr('network')
@attr('slow')
def test_quick_model():
    uniprot_entry_name = 'EGFR_HUMAN'
    uniprot_domain_regex = '^Protein kinase'
    pdbids = ['4AF3']
    chainids = {
        '4AF3': ['A']
    }
    with enter_temp_dir():
        QuickModel(
            target_uniprot_entry_name=uniprot_entry_name,
            uniprot_domain_regex=uniprot_domain_regex,
            pdbids=pdbids,
            chainids=chainids,
            loopmodel=False,
            sim_length=2.0*unit.femtoseconds,
        )
        assert os.path.exists(os.path.join(
            default_project_dirnames.models,
            'EGFR_HUMAN_D0', 'AURKB_HUMAN_D0_4AF3_A', 'explicit-refined.pdb.gz'
        ))
