import os
import tempfile
import shutil
import distutils.dir_util
import contextlib
import ensembler.initproject


@contextlib.contextmanager
def integration_test_context(set_up_project_stage='init'):
    temp_dir = tempfile.mkdtemp()
    ensembler.initproject.InitProject(temp_dir)
    cwd = os.getcwd()

    set_up_project_stage_methods = SetUpProjectStageMethods(temp_dir)
    set_up_method = getattr(set_up_project_stage_methods, set_up_project_stage)
    set_up_method()

    os.chdir(temp_dir)
    yield temp_dir
    os.chdir(cwd)
    shutil.rmtree(temp_dir)


class SetUpProjectStageMethods:
    def __init__(self, temp_dir):
        self.temp_dir = temp_dir

    def init(self):
        shutil.copy(os.path.join('tests', 'integration_test_resources', 'meta0.yaml'), self.temp_dir)

    def targets(self):
        shutil.copy(os.path.join('tests', 'integration_test_resources', 'meta0.yaml'), self.temp_dir)
        distutils.dir_util.copy_tree(os.path.join('tests', 'integration_test_resources', 'targets'), os.path.join(self.temp_dir, 'targets'))

    def templates_resolved(self):
        shutil.copy(os.path.join('tests', 'integration_test_resources', 'meta0.yaml'), self.temp_dir)
        distutils.dir_util.copy_tree(os.path.join('tests', 'integration_test_resources', 'targets'), os.path.join(self.temp_dir, 'targets'))
        shutil.copy(os.path.join('tests', 'integration_test_resources', 'templates', 'meta0.yaml'), os.path.join(self.temp_dir, 'templates'))
        shutil.copy(os.path.join('tests', 'integration_test_resources', 'templates', 'templates-resolved.fa'), os.path.join(self.temp_dir, 'templates'))
        shutil.copy(os.path.join('tests', 'integration_test_resources', 'templates', 'templates-full-seq.fa'), os.path.join(self.temp_dir, 'templates'))
        distutils.dir_util.copy_tree(os.path.join('tests', 'integration_test_resources', 'templates', 'structures-resolved'), os.path.join(self.temp_dir, 'templates', 'structures-resolved'))

    def templates_modeled_loops(self):
        shutil.copy(os.path.join('tests', 'integration_test_resources', 'meta0.yaml'), self.temp_dir)
        distutils.dir_util.copy_tree(os.path.join('tests', 'integration_test_resources', 'targets'), os.path.join(self.temp_dir, 'targets'))
        shutil.copy(os.path.join('tests', 'integration_test_resources', 'templates', 'meta0.yaml'), os.path.join(self.temp_dir, 'templates'))
        shutil.copy(os.path.join('tests', 'integration_test_resources', 'templates', 'templates-resolved.fa'), os.path.join(self.temp_dir, 'templates'))
        shutil.copy(os.path.join('tests', 'integration_test_resources', 'templates', 'templates-full-seq.fa'), os.path.join(self.temp_dir, 'templates'))
        distutils.dir_util.copy_tree(os.path.join('tests', 'integration_test_resources', 'templates', 'structures-resolved'), os.path.join(self.temp_dir, 'templates', 'structures-resolved'))
        distutils.dir_util.copy_tree(os.path.join('tests', 'integration_test_resources', 'templates', 'structures-modeled-loops'), os.path.join(self.temp_dir, 'templates', 'structures-modeled-loops'))