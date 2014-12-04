import os
import tempfile
import shutil
import distutils.dir_util
import contextlib
import msmseeder.initproject


@contextlib.contextmanager
def integration_test_context(project_data='init'):
    temp_dir = tempfile.mkdtemp()
    msmseeder.initproject.initproject(temp_dir)
    cwd = os.getcwd()

    shutil.copy('meta0.yaml', temp_dir)
    if project_data == 'templates_resolved':
        shutil.copy('meta0.yaml', temp_dir)
        distutils.dir_util.copy_tree(os.path.join('targets'), os.path.join(temp_dir, 'targets'))
        shutil.copy(os.path.join('templates', 'meta0.yaml'), os.path.join(temp_dir, 'templates'))
        shutil.copy(os.path.join('templates', 'templates-resolved.fa'), os.path.join(temp_dir, 'templates'))
        shutil.copy(os.path.join('templates', 'templates-full-seq.fa'), os.path.join(temp_dir, 'templates'))
        distutils.dir_util.copy_tree(os.path.join('templates', 'structures-resolved'), os.path.join(temp_dir, 'templates', 'structures-resolved'))

    os.chdir(temp_dir)
    yield temp_dir
    os.chdir(cwd)
    shutil.rmtree(temp_dir)