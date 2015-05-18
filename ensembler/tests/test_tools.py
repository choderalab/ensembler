from ensembler.tools.inspect import LoopmodelLogs
from ensembler.tests.integrationtest_utils import integrationtest_context


def test_loopmodel_logs():
    with integrationtest_context(set_up_project_stage='templates_modeled_loops'):
        loopmodel_logs = LoopmodelLogs()
        loopmodel_logs.add_missing_resis_data()