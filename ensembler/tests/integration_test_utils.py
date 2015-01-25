import os
import tempfile
import shutil
import distutils.dir_util
import contextlib
import ensembler.initproject

integration_test_resources_dir = os.path.join('tests', 'integration_test_resources')


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
        self.targets_list = ['EGFR_HUMAN_D0', 'KC1D_HUMAN_D0']
        self.templates_list = ['KC1D_HUMAN_D0_4KB8_D', 'KC1D_HUMAN_D0_4HNF_A']

    def init(self):
        shutil.copy(os.path.join(integration_test_resources_dir, 'meta0.yaml'), self.temp_dir)

    def targets(self):
        self.init()
        distutils.dir_util.copy_tree(os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.targets), os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.targets))

    def templates_resolved(self):
        self.targets()
        shutil.copy(os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.templates, 'meta0.yaml'), os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.templates))
        shutil.copy(os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.templates, 'templates-resolved-seq.fa'), os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.templates))
        shutil.copy(os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.templates, 'templates-full-seq.fa'), os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.templates))
        distutils.dir_util.copy_tree(os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.templates_structures_resolved), os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.templates_structures_resolved))

    def templates_modeled_loops(self):
        self.templates_resolved()
        distutils.dir_util.copy_tree(os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.templates_structures_modeled_loops), os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.templates_structures_modeled_loops))

    def aligned(self):
        self.templates_modeled_loops()
        for target in self.targets_list:
            for template in self.templates_list:
                ensembler.utils.create_dir(os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.models, target, template))
                shutil.copy(
                    os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.models, target, template, 'alignment.pir'),
                    os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.models, target, template, 'alignment.pir')
                )
            shutil.copy(
                os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.models, target, 'sequence-identities.txt'),
                os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.models, target)
            )

    def modeled(self):
        self.aligned()
        for target in self.targets_list:
            shutil.copy(
                os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.models, target, 'build_models-meta0.yaml'),
                os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.models, target)
            )
            for template in self.templates_list:
                for filename in [
                    'model.pdb',
                    'model.pdb.gz',
                    'modeling-log.yaml',
                    'restraints.rsr.gz',
                    'sequence-identity.txt',
                ]:
                    shutil.copy(
                        os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.models, target, template, filename),
                        os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.models, target, template)
                    )

    def clustered(self):
        self.modeled()
        for target in self.targets_list:
            for filename in [
                'cluster_models-meta0.yaml',
                'unique-models.txt',
            ]:
                shutil.copy(
                    os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.models, target, filename),
                    os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.models, target)
                )
            for template in self.templates_list:
                shutil.copy(
                    os.path.join(integration_test_resources_dir, ensembler.core.default_project_dirnames.models, target, template, 'unique_by_clustering'),
                    os.path.join(self.temp_dir, ensembler.core.default_project_dirnames.models, target, template)
                )