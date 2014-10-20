import os
import shutil
import msmseeder
import msmseeder.tests
import msmseeder.modelling
from mock import Mock


def test_build_model():
    template_filepath = os.path.abspath(os.path.join('tests', 'resources', 'mock_template.pdb'))

    with msmseeder.tests.utils.enter_temp_directory():
        # given
        target = Mock()
        template = Mock()
        target.id = 'mock_target'
        target.seq = 'YILGDTLGVGGKVKVGKH'
        template.id = 'mock_template'
        template.seq = 'YQNLSPVGSGAYGSVCAAFD'

        shutil.copy(template_filepath, '.')

        # when
        output_text = msmseeder.modelling.build_model(
            target,
            template,
            template_structure_dir='.',
            aln_filepath='alignment.pir',
            seqid_filepath='sequence-identity.txt',
            model_pdbfilepath='model.pdb.gz',
            restraint_filepath='restraints.rsr.gz'
        )

        # then

        # Example model.pdb.gz contents (not testing this as it may be dependent
        # upon Modeller version as well as modelling stochasticity):
        #
        # ..EXPDTA    THEORETICAL MODEL, MODELLER 9.12 2014/08/26 13:15:44
        # REMARK   6 MODELLER OBJECTIVE FUNCTION:       326.6798
        # REMARK   6 MODELLER BEST TEMPLATE % SEQ ID:  27.778
        # ATOM      1  N   TYR     1      48.812  50.583  13.949  1.00110.28           N
        # ATOM      2  CA  TYR     1      49.070  50.334  15.387  1.00110.28           C

        assert os.path.exists('model.pdb.gz')
        assert os.path.getsize('model.pdb.gz') > 0