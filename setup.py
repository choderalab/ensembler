import os
import subprocess
from setuptools import setup

##########################
VERSION = "0.2"
ISRELEASED = False
__version__ = VERSION
##########################

def read_readme(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

##########################
# Function for determining current git commit
##########################

def git_version():
    # Return the git revision as a string
    # copied from numpy setup.py
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = 'Unknown'

    return GIT_REVISION

##########################
# Function for writing version.py (this will be copied to the install directory)
##########################

msmseeder_version_filepath = 'msmseeder/version.py'
def write_version_py(filename=msmseeder_version_filepath):
    cnt = """# THIS FILE IS GENERATED FROM MSMSEEDER SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    # Adding the git rev number needs to be done inside write_version_py(),
    # otherwise the import of numpy.version messes up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    else:
        GIT_REVISION = 'Unknown'

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()

##########################
# Find packages
##########################

def find_packages():
    """Find all of msmseeder's python packages.
    Adapted from IPython's setupbase.py. Copyright IPython
    contributors, licensed under the BSD license.
    """
    packages = ['msmseeder.scripts']
    for dir,subdirs,files in os.walk('msmseeder'):
        package = dir.replace(os.path.sep, '.')
        if '__init__.py' not in files:
            # not a package
            continue
        packages.append(package)
    return packages

##########################
# Setup
##########################

write_version_py()
setup(
    name = 'msmseeder',
    version = __version__,
    author = 'Daniel L Parton',
    author_email = 'daniel.parton@choderalab.org',
    description = 'Generation of diverse protein structural ensembles, for the initialization of molecular dynamics simulations and subsequent construction of Markov state models. ',
    license='GPLv2',
    long_description = read_readme('README.md'),
    packages = find_packages(),
    package_dir = {'msmseeder.scripts': 'scripts'},
    package_data = {'msmseeder.tests': ['resources/*']},
    # scripts = ['scripts/MSMSeederInit.py', 'scripts/MSMSeederGatherTargets.py', 'scripts/MSMSeederGatherTemplates.py', 'scripts/MSMSeederBuildModels.py', 'scripts/MSMSeederRefineImplicitMD.py', 'scripts/MSMSeederSolvate.py', 'scripts/MSMSeederRefineExplicitMD.py', 'scripts/MSMSeederPackageModels.py'],
    # data_files = [('', ['LICENSE']), ('templates', ['project-data.yaml-TEMPLATE', 'manual-specifications.yaml-TEMPLATE'])],
    entry_points = {'console_scripts':
        ['MSMSeederInit = msmseeder.scripts.MSMSeederInit:main']
    },
)

# Delete version.py file from this directory
os.remove(msmseeder_version_filepath)
