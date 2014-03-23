import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'MSMSeeder',
    version = '0.2',
    author = 'Daniel L Parton',
    author_email = 'daniel.parton@choderalab.org',
    description = 'Generation of diverse protein structural ensembles, for the initialization of molecular dynamics simulations and subsequent construction of Markov state models. ',
    long_description = read('README.md'),
    packages = ['MSMSeeder', 'tests'],
    scripts = ['scripts/InitMSMSeederProject.py', 'scripts/GatherTargets.py', 'scripts/GatherTemplates.py', 'scripts/BuildModels.py', 'scripts/RefineImplicitMD.py', 'scripts/Solvate.py', 'scripts/RefineExplicitMD.py'],
    data_files = [('', ['LICENSE']), ('templates', ['project-data.yaml-TEMPLATE', 'manual-specifications.yaml-TEMPLATE'])],
)
