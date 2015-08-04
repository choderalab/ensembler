import os
import tempfile
import shutil
import distutils.dir_util
import contextlib
import ensembler.initproject
from ensembler.core import default_project_dirnames
from ensembler.tests.utils import get_installed_resource_filename


@contextlib.contextmanager
def integrationtest_context(set_up_project_stage='init'):
    temp_dir = tempfile.mkdtemp()
    ensembler.initproject.InitProject(temp_dir)
    cwd = os.getcwd()

    set_up_project_stage_methods = SetUpSampleProject(temp_dir)
    set_up_method = getattr(set_up_project_stage_methods, set_up_project_stage)
    set_up_method()

    os.chdir(temp_dir)
    yield temp_dir
    os.chdir(cwd)
    shutil.rmtree(temp_dir)


class SetUpSampleProject:
    def __init__(self, project_dir):
        self.project_dir = project_dir
        self.targets_list = ['EGFR_HUMAN_D0', 'KC1D_HUMAN_D0']
        self.templates_list = ['KC1D_HUMAN_D0_4KB8_D', 'KC1D_HUMAN_D0_4HNF_A']

    def init(self):
        ensembler.initproject.InitProject(self.project_dir)
        shutil.copy(get_installed_resource_filename(os.path.join('example_project', 'meta0.yaml')), self.project_dir)

    def targets(self):
        self.init()
        distutils.dir_util.copy_tree(
            get_installed_resource_filename(
                os.path.join('example_project', default_project_dirnames.targets)
            ),
            os.path.join(self.project_dir, default_project_dirnames.targets)
        )

    def templates_resolved(self):
        self.targets()
        shutil.copy(
            get_installed_resource_filename(
                os.path.join('example_project', default_project_dirnames.templates, 'meta0.yaml')
            ),
            os.path.join(self.project_dir, default_project_dirnames.templates))
        shutil.copy(
            get_installed_resource_filename(
                os.path.join('example_project', default_project_dirnames.templates, 'templates-resolved-seq.fa')
            ),
            os.path.join(self.project_dir, default_project_dirnames.templates))
        shutil.copy(
            get_installed_resource_filename(
                os.path.join('example_project', default_project_dirnames.templates, 'templates-full-seq.fa')
            ),
            os.path.join(self.project_dir, default_project_dirnames.templates)
        )
        distutils.dir_util.copy_tree(
            get_installed_resource_filename(
                os.path.join('example_project', default_project_dirnames.templates_structures_resolved)
            ),
            os.path.join(self.project_dir, default_project_dirnames.templates_structures_resolved)
        )

    def templates_modeled_loops(self):
        self.templates_resolved()
        distutils.dir_util.copy_tree(
            get_installed_resource_filename(
                os.path.join('example_project', default_project_dirnames.templates_structures_modeled_loops)
            ),
            os.path.join(self.project_dir, default_project_dirnames.templates_structures_modeled_loops)
        )

    def aligned(self):
        self.templates_modeled_loops()
        for target in self.targets_list:
            for template in self.templates_list:
                ensembler.utils.create_dir(
                    os.path.join(self.project_dir, default_project_dirnames.models, target, template)
                )
        self._copy_modeling_files(
            target_level_files=[
                'sequence-identities.txt'
            ],
            template_level_files=[
                'alignment.pir'
            ]
        )

    def modeled(self):
        self.aligned()
        self._copy_modeling_files(
            target_level_files=[
                'build_models-meta0.yaml'
            ],
            template_level_files=[
                'model.pdb',
                'model.pdb.gz',
                'modeling-log.yaml',
                'restraints.rsr.gz',
                'sequence-identity.txt',
            ]
        )

    def clustered(self):
        self.modeled()
        self._copy_modeling_files(
            target_level_files=[
                'cluster_models-meta0.yaml',
                'unique-models.txt',
            ],
            template_level_files=[
                'unique_by_clustering'
            ]
        )

    def refined_implicit(self):
        self.clustered()
        self._copy_modeling_files(
            target_level_files=[
                'refine_implicit_md-meta0.yaml'
            ],
            template_level_files=[
                'implicit-refined.pdb.gz',
                'implicit-log.yaml',
                'implicit-energies.txt',
            ]
        )

    def solvated(self):
        self.refined_implicit()
        self._copy_modeling_files(
            target_level_files=[
                'solvate_models-meta0.yaml',
                'determine_nwaters-meta0.yaml',
                'nwaters-use.txt',
                'nwaters-max.txt',
                'nwaters-list.txt',
            ],
            template_level_files=[
                'nwaters.txt',
            ]
        )

    def refined_explicit(self):
        self.solvated()
        self._copy_modeling_files(
            target_level_files=[
                'refine_explicit_md-meta0.yaml',
            ],
            template_level_files=[
                'explicit-system.xml.gz',
                'explicit-state.xml.gz',
                'explicit-refined.pdb.gz',
                'explicit-log.yaml',
                'explicit-integrator.xml.gz',
                'explicit-energies.txt',
            ]
        )

    def _copy_modeling_files(self, target_level_files=None, template_level_files=None):
        for target in self.targets_list:
            for filename in target_level_files:
                shutil.copy(
                    get_installed_resource_filename(os.path.join('example_project', default_project_dirnames.models, target, filename)),
                    os.path.join(self.project_dir, default_project_dirnames.models, target)
                )
            for template in self.templates_list:
                for filename in template_level_files:
                    shutil.copy(
                        get_installed_resource_filename(os.path.join('example_project', default_project_dirnames.models, target, template, filename)),
                        os.path.join(self.project_dir, default_project_dirnames.models, target, template)
                    )