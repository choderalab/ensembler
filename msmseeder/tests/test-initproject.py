import os
import shutil
import tempfile
import msmseeder.initproject

def test_initproject():
    tmpdir = tempfile.mkdtemp()
    msmseeder.initproject.initproject(tmpdir)
    for dirpath in ['targets', 'structures', 'templates', 'models', 'packaged-models',
                        'structures/pdb',
                        'structures/sifts',
                        'templates/structures']:
        assert os.path.exists(os.path.join(tmpdir, dirpath))
    assert os.path.exists(os.path.join(tmpdir, 'meta.yaml'))
    shutil.rmtree(tmpdir)

# TODO def test_gather_targets_from_targetexplorer():

# TODO def test_gather_targets_from_uniprot():