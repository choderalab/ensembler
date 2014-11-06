import os
import shutil
import msmseeder
import msmseeder.tests
import msmseeder.modelling
from mock import Mock
import datetime
from msmseeder.utils import enter_temp_dir


def test_build_model():
    template_filepath = os.path.abspath(os.path.join('tests', 'resources', 'mock_template.pdb'))

    with enter_temp_dir():
        target = Mock()
        template = Mock()
        target_setup_data = Mock()
        target.id = 'mock_target'
        target.seq = 'YILGDTLGVGGKVKVGKH'
        template.id = 'mock_template'
        template.seq = 'YQNLSPVGSGAYGSVCAAFD'
        target_setup_data.target_starttime = datetime.datetime.utcnow()
        target_setup_data.models_target_dir = os.path.join(msmseeder.core.default_project_dirnames.models, target.id)
        os.mkdir(msmseeder.core.default_project_dirnames.models)
        os.mkdir(msmseeder.core.default_project_dirnames.templates)
        os.mkdir(msmseeder.core.default_project_dirnames.templates_structures)
        os.mkdir(target_setup_data.models_target_dir)

        shutil.copy(template_filepath, msmseeder.core.default_project_dirnames.templates_structures)

        msmseeder.modelling.build_model(
            target,
            template,
            target_setup_data=target_setup_data
        )

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