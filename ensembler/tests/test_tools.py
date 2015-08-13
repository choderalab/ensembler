from ensembler.tools.inspect import LoopmodelLogs
from ensembler.tools.mktraj import MkTrajImplicitStart
from ensembler.tests.integrationtest_utils import integrationtest_context


def test_loopmodel_logs():
    with integrationtest_context(set_up_project_stage='templates_modeled_loops'):
        loopmodel_logs = LoopmodelLogs()
        loopmodel_logs.add_missing_resis_data()


def test_mktraj_implicit_start():
    with integrationtest_context(set_up_project_stage='refined_explicit'):
        MkTrajImplicitStart(targetid='EGFR_HUMAN_D0', loglevel='debug')
