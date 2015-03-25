import os
import logging
import sys
import re
import warnings
import numpy as np
import Bio
import Bio.SeqIO
from collections import namedtuple

# ========
# Global package variables
# ========
import ensembler

installation_toplevel_dir = os.path.abspath(os.path.join(os.path.dirname(__file__)))

datestamp_format_string = '%Y-%m-%d %H:%M:%S UTC'

manual_overrides_filename = 'manual-overrides.yaml'

template_acceptable_ratio_resolved_residues = 0.7

# listed in order
project_stages = [
    'init',
    'gather_targets',
    'gather_templates',
    'build_models',
    'cluster_models',
    'refine_implicit_md',
    'solvate_models',
    'determine_nwaters',
    'refine_explicit_md',
    'package_for_fah',
]

project_dirtypes = [
    'targets',
    'templates',
    'structures', 'models',
    'packaged_models',
    'structures_pdb',
    'structures_sifts',
    'templates_structures_resolved',
    'templates_structures_modeled_loops',
]
ProjectDirNames = namedtuple('ProjectDirNames', project_dirtypes)
default_project_dirnames = ProjectDirNames(
    targets='targets',
    templates='templates',
    structures='structures',
    models='models',
    packaged_models='packaged_models',
    structures_pdb=os.path.join('structures', 'pdb'),
    structures_sifts=os.path.join('structures', 'sifts'),
    templates_structures_resolved=os.path.join('templates', 'structures-resolved'),
    templates_structures_modeled_loops=os.path.join('templates', 'structures-modeled-loops'),
)

installation_dirtypes = [
    'tests_integration_test_resources'
]
InstallationDirNames = namedtuple('InstallationDirNames', installation_dirtypes)
default_installation_dirnames = InstallationDirNames(
    tests_integration_test_resources=os.path.join('tests', 'integration_test_resources'),
)

logger = logging.getLogger('info')
default_loglevel = 'info'
loglevel_obj = getattr(logging, default_loglevel.upper())
logger.setLevel(loglevel_obj)
logger.addHandler(logging.StreamHandler(stream=sys.stdout))

# metadata_filename_regex = re.compile('(^|.*)(meta)([0-9]+)\.yaml')

target_id_regex = re.compile('^[A-Z0-9]{2,5}_[A-Z0-9]{2,5}_D[0-9]{1,3}$')
template_id_regex = re.compile('^[A-Z0-9]{2,5}_[A-Z0-9]{2,5}_D[0-9]{1,3}_[A-Z0-9]{4}_[A-Z0-9]$')

model_filenames_by_ensembler_stage = {
    'build_models': 'model.pdb.gz',
    'refine_implicit_md': 'implicit-refined.pdb.gz',
    'refine_explicit_md': 'explicit-refined.pdb.gz',
}

# ========
# MPI
# ========

class MPIState:
    def __init__(self):
        try:
            import mpi4py.MPI
            self.comm = mpi4py.MPI.COMM_WORLD
            self.rank = self.comm.rank
            self.size = self.comm.size
        except Exception as e:
            logger.debug('Error initializing MPIState:\n%s' % e)
            self.comm = None
            self.rank = 0
            self.size = 1

mpistate = MPIState()

# ========
# YAML
# ========

try:
    from yaml import CLoader as YamlLoader, CDumper as YamlDumper
except ImportError:
    from yaml import Loader as YamlLoader, Dumper as YamlDumper

class literal_str(str): pass

def change_style(style, representer):
    def new_representer(dumper, data):
        scalar = representer(dumper, data)
        scalar.style = style
        return scalar
    return new_representer

from yaml.representer import SafeRepresenter

# represent_str does handle some corner cases, so use that
# instead of calling represent_scalar directly
represent_literal_str = change_style('|', SafeRepresenter.represent_str)

import yaml
yaml.add_representer(literal_str, represent_literal_str)

# ========
# Definitions
# ========

def get_utcnow_formatted():
    import datetime
    now = datetime.datetime.utcnow()
    datestamp = now.strftime(datestamp_format_string)
    return datestamp

def strf_timedelta(delta_t):
    hours, remainder = divmod(delta_t.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return '%d:%d:%d' % (hours, minutes, seconds)

def check_project_toplevel_dir(raise_exception=True):
    import os
    for dirtype in project_dirtypes:
        if dirtype == 'packaged_models':
            continue
        dirpath = getattr(default_project_dirnames, dirtype)
        if not os.path.exists(dirpath):
            warnings.warn('Directory %s not found' % dirpath)
            if raise_exception:
                raise Exception('This is not the top-level directory of an Ensembler project.')
            else:
                return False
    return True

class LogFile:
    def __init__(self, log_filepath):
        import socket
        self.log_filepath = log_filepath
        self.log_data = {
            'datestamp': get_utcnow_formatted(),
            'hostname': socket.gethostname(),
        }

    def log(self, new_log_data={}):
        self.log_data.update(new_log_data)

        with open(self.log_filepath, 'w') as log_file:
            yaml.dump(self.log_data, log_file, default_flow_style=False, Dumper=YamlDumper)


class ManualOverrides:
    def __init__(self):
        import yaml
        if os.path.exists(manual_overrides_filename):
            with open(manual_overrides_filename, 'r') as manual_overrides_file:
                manual_overrides_yaml = yaml.load(manual_overrides_file, Loader=YamlLoader)
        else:
            manual_overrides_yaml = {}

        self.target = TargetManualOverrides(manual_overrides_yaml)
        self.template = TemplateManualOverrides(manual_overrides_yaml)


class TargetManualOverrides:
    '''
    Parameters
    ----------
    manual_overrides_yaml: dict

    Attributes
    ----------
    domain_spans: dict
        dict with structure {`targetid`: `domain_span`, ...} where
        `domain_span` is a str e.g. '242-513'
    '''
    def __init__(self, manual_overrides_yaml):
        target_dict = manual_overrides_yaml.get('target-selection')
        if target_dict != None:
            self.domain_spans = target_dict.get('domain-spans')
        else:
            self.domain_spans = {}


class TemplateManualOverrides:
    '''
    Parameters
    ----------
    manual_overrides_yaml: dict

    Attributes
    ----------
    min_domain_len: int or NoneType
    max_domain_len: int or NoneType
    domain_spans: dict
        dict with structure {`targetid`: `domain_span`, ...} where
        `domain_span` is a str e.g. '242-513'
    skip_pdbs: list
        list of PDB IDs to skip
    '''
    def __init__(self, manual_overrides_yaml):
        template_dict = manual_overrides_yaml.get('template-selection')
        if template_dict != None:
            self.min_domain_len = template_dict.get('min-domain-len')
            self.max_domain_len = template_dict.get('max-domain-len')
            self.domain_spans = template_dict.get('domain-spans')
            self.skip_pdbs = template_dict.get('skip-pdbs')
        else:
            self.min_domain_len = None
            self.max_domain_len = None
            self.domain_spans = {}
            self.skip_pdbs = []


def gen_metadata_filename(ensembler_stage, metadata_file_index):
    for modelling_stage in ['build_models', 'cluster_models', 'refine_implicit_md', 'solvate_models', 'determine_nwaters', 'refine_explicit_md']:
        if ensembler_stage == modelling_stage:
            return '%s-meta%d.yaml' % (ensembler_stage, metadata_file_index)
    return 'meta%d.yaml' % metadata_file_index


class ProjectMetadata:
    # TODO - have metadata files output as follows
    # models/SRC_HUMAN_D0/build_models-meta.yaml
    # models/SRC_HUMAN_D0/implicit_refinement-meta.yaml
    """
    Examples
    --------
    >>> init_project_metadata = ensembler.core.ProjectMetadata(project_stage='init')
    >>> test_data = {'test_field': 'test_value'}
    >>> init_project_metadata.add_data(test_data)
    >>> init_project_metadata.write()
    >>> gather_targets_project_metadata = ensembler.core.ProjectMetadata(project_stage='gather_targets')
    >>> gather_targets_project_metadata.add_data(test_data)
    >>> gather_targets_project_metadata.write()
    """
    def __init__(self, project_stage='init', target_id=None, project_toplevel_dir='.'):
        self.data = {}
        self.project_stage = project_stage
        self.target_id = target_id
        self.project_toplevel_dir = project_toplevel_dir
        if project_stage != project_stages[0]:
            self.add_all_prev_metadata(project_stage)

    def add_all_prev_metadata(self, project_stage):
        current_project_stage_index = project_stages.index(project_stage)
        for prev_project_stage in project_stages[: current_project_stage_index]:
            self.add_prev_metadata(prev_project_stage)

    def add_prev_metadata(self, project_stage):
        latest_metadata_filepath = self.determine_latest_metadata_filepath(project_stage)
        with open(latest_metadata_filepath) as latest_metadata_file:
            prev_metadata = yaml.load(latest_metadata_file, Loader=YamlLoader)
        self.add_data(prev_metadata[project_stage], project_stage=project_stage)

    def determine_latest_metadata_filepath(self, project_stage):
        metadata_dir = self.metadata_dir_mapper(project_stage, self.target_id)
        metadata_file_basename = self.metadata_file_basename_mapper(project_stage)
        latest_metadata_file_index = self.determine_latest_metadata_file_index(project_stage)
        latest_metadata_filepath = os.path.join(metadata_dir, '%s%d.yaml' % (metadata_file_basename, latest_metadata_file_index))
        return latest_metadata_filepath

    def metadata_dir_mapper(self, project_stage, target_id=None):
        metadata_dir_dict = {
            'init': '.',
            'gather_targets': 'targets',
            'gather_templates': 'templates',
        }
        if project_stage in metadata_dir_dict:
            return metadata_dir_dict[project_stage]
        elif project_stage in ['build_models', 'cluster_models', 'refine_implicit_md', 'solvate_models', 'determine_nwaters', 'refine_explicit_md']:
            return os.path.join('models', target_id)

    def metadata_file_basename_mapper(self, project_stage):
        if project_stage in ['build_models', 'cluster_models', 'refine_implicit_md', 'solvate_models', 'determine_nwaters', 'refine_explicit_md']:
            return '%s-meta' % project_stage
        else:
            return 'meta'

    def determine_latest_metadata_file_index(self, project_stage):
        """
        Returns -1 if no metadata files found
        :param project_stage: str
        :return: int
        """
        metadata_dir = self.metadata_dir_mapper(project_stage, target_id=self.target_id)
        dir_contents = os.listdir(metadata_dir)
        metadata_file_basename = self.metadata_file_basename_mapper(project_stage)
        metadata_file_regex = re.compile('%s([0-9]+)\.yaml' % metadata_file_basename)
        metadata_file_indices = []
        for filename in dir_contents:
            match = re.search(metadata_file_regex, filename)
            if match:
                metadata_file_indices.append(int(match.groups()[0]))
        if len(metadata_file_indices) > 0:
            return max(metadata_file_indices)
        else:
            return -1

    def gen_metadata_filepath_from_dir_index_and_file_basename(self, dirpath, file_basename, index):
        metadata_filepath = os.path.join(dirpath, '%s%d.yaml' % (file_basename, index))
        return metadata_filepath

    def add_data(self, data, project_stage=None):
        """
        Add metadata to the ProjectMetadata object.
        If project_stage is not passed as an argument, it is set to ProjectMetadata.project_stage (which is itself set during object initialization).

        Parameters
        ----------
        data: dict
        project_stage: str

        Examples
        --------
        >>> project_metadata = ensembler.core.ProjectMetadata(project_stage='init')
        >>> metadata = {'datestamp': ensembler.core.get_utcnow_formatted()}
        >>> project_metadata.add_data(metadata)
        """
        if project_stage is None:
            project_stage = self.project_stage
        self.data.update({
            project_stage: data
        })

    def add_iteration_number_to_metadata(self, iter_number):
        """
        Parameters
        ----------
        iter_number: int
        """
        data = self.data[self.project_stage]
        data['iteration'] = iter_number

    def write(self):
        metadata_dir = os.path.join(self.project_toplevel_dir, self.metadata_dir_mapper(self.project_stage, target_id=self.target_id))
        metadata_file_basename = self.metadata_file_basename_mapper(self.project_stage)
        latest_metadata_file_index = self.determine_latest_metadata_file_index(self.project_stage)
        self.add_iteration_number_to_metadata(latest_metadata_file_index+1)
        metadata_filepath = self.gen_metadata_filepath_from_dir_index_and_file_basename(metadata_dir, metadata_file_basename, latest_metadata_file_index+1)
        with open(metadata_filepath, 'w') as ofile:
            for stage in project_stages:
                if stage in self.data.keys():
                    subdict = {stage: self.data[stage]}
                    yaml.dump(subdict, ofile, default_flow_style=False, Dumper=YamlDumper)


class DeprecatedProjectMetadata:
    # TODO deprecate
    def __init__(self, data):
        self.data = data

    def write(self, ofilepath):
        with open(ofilepath, 'w') as ofile:
            for stage in project_stages:
                if stage in self.data.keys():
                    subdict = {stage: self.data[stage]}
                    yaml.dump(subdict, ofile, default_flow_style=False, Dumper=YamlDumper)


def write_metadata(new_metadata_dict, ensembler_stage, target_id=None):
    # TODO deprecate
    if ensembler_stage == 'init':
        metadata_dict = {}
    else:
        prev_ensembler_stage = project_stages[project_stages.index(ensembler_stage) - 1]
        prev_metadata_filepath = metadata_file_mapper(prev_ensembler_stage, target_id=target_id)
        with open(prev_metadata_filepath) as prev_metadata_file:
            metadata_dict = yaml.load(prev_metadata_file, Loader=YamlLoader)

    metadata_dict.update(new_metadata_dict)
    metadata = ProjectMetadata(metadata_dict)
    metadata.write(metadata_file_mapper(ensembler_stage, target_id=target_id))


def metadata_file_mapper(ensembler_stage, target_id=None):
    # TODO deprecate
    metadata_file_dict = {
        'init': 'meta.yaml',
        'gather_targets': os.path.join('targets', 'meta.yaml'),
        'gather_templates': os.path.join('templates', 'meta.yaml'),
    }
    if ensembler_stage in metadata_file_dict:
        return metadata_file_dict[ensembler_stage]
    elif ensembler_stage in ['build_models', 'cluster_models', 'refine_implicit_md', 'solvate_models', 'determine_nwaters', 'refine_explicit_md']:
        return os.path.join('models', target_id, 'meta.yaml')


def encode_url_query(uniprot_query):
    def replace_all(text, replace_dict):
        for i, j in replace_dict.iteritems():
            text = text.replace(i, j)
        return text

    encoding_dict = {
        ' ': '+',
        ':': '%3A',
        '(': '%28',
        ')': '%29',
        '"': '%22',
        '=': '%3D',
    }
    return replace_all(uniprot_query, encoding_dict)

def xpath_match_regex_case_sensitive(context, attrib_values, xpath_argument):
    ''' To be used as an lxml XPath extension, for regex searches of attrib values.
    '''
    import re
    # If no attrib found
    if len(attrib_values) == 0:
        return False
    # If attrib found, then run match against regex
    else:
        regex = re.compile(xpath_argument)
        return bool( re.search(regex, attrib_values[0]) )

def xpath_match_regex_case_insensitive(context, attrib_values, xpath_argument):
    ''' To be used as an lxml XPath extension, for regex searches of attrib values.
    '''
    import re
    # If no attrib found
    if len(attrib_values) == 0:
        return False
    # If attrib found, then run match against regex
    else:
        regex = re.compile(xpath_argument, re.IGNORECASE)
        return bool( re.search(regex, attrib_values[0]) )

def sequnwrap(sequence):
    '''
    Unwraps a wrapped sequence string
    '''
    unwrapped = sequence.strip()
    unwrapped = ''.join(unwrapped.split('\n'))
    return unwrapped

def seqwrap(sequence, add_star=False):
    '''
    Wraps a sequence string to a width of 60.
    If add_star is set to true, an asterisk will be added
    to the end of the sequence, for compatibility with
    Modeller.
    '''
    if add_star:
        sequence += '*'
    wrapped = ''
    for i in range(0,len(sequence),60):
        wrapped += sequence[i: i+60] + '\n'
    return wrapped


def construct_fasta_str(id, seq):
    target_fasta_string = '>%s\n%s\n' % (id, seqwrap(seq).strip())
    return target_fasta_string


def get_targets():
    targets_dir = os.path.abspath('targets')
    targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
    targets = list(Bio.SeqIO.parse(targets_fasta_filename, 'fasta'))
    return targets


def get_templates():
    templates_dir = os.path.abspath('templates')
    templates_resolved_seq_fasta_filename = os.path.join(templates_dir, 'templates-resolved-seq.fa')
    templates_resolved_seq = list(Bio.SeqIO.parse(templates_resolved_seq_fasta_filename, 'fasta'))
    templates_full_seq_fasta_filename = os.path.join(templates_dir, 'templates-full-seq.fa')
    templates_full_seq = list(Bio.SeqIO.parse(templates_full_seq_fasta_filename, 'fasta'))
    return templates_resolved_seq, templates_full_seq


def get_targets_and_templates():
    targets = get_targets()
    templates_resolved_seq, templates_full_seq = get_templates()
    return targets, templates_resolved_seq, templates_full_seq


def find_loopmodel_executable():
    for path in os.environ['PATH'].split(os.pathsep):
        if not os.path.exists(path):
            continue
        path = path.strip('"')
        for filename in os.listdir(path):
            if len(filename) >= 10 and filename[0: 10] == 'loopmodel.':
                if filename[-5:] == 'debug':
                    warnings.warn('loopmodel debug version (%s) will be ignored - runs extremely slow' % filename)
                    continue
                return os.path.join(path, filename)
    raise Exception('Loopmodel executable not found in PATH')


def check_ensembler_modeling_stage_first_model_file_exists(ensembler_stage, targetid):
    models_target_dir = os.path.join(ensembler.core.default_project_dirnames.models, targetid)
    root, dirnames, filenames = next(os.walk(models_target_dir))
    for dirname in dirnames:
        if re.match(ensembler.core.template_id_regex, dirname):
            model_filepath = os.path.join(
                models_target_dir,
                dirname,
                ensembler.core.model_filenames_by_ensembler_stage[ensembler_stage]
            )
            if os.path.exists(model_filepath):
                return True
    return False


def check_ensembler_modeling_stage_metadata_exists(ensembler_stage, targetid):
    next_ensembler_stage = ensembler.core.project_stages[
        ensembler.core.project_stages.index(ensembler_stage) + 1]
    try:
        project_metadata = ensembler.core.ProjectMetadata(project_stage=next_ensembler_stage, target_id=targetid)
    except IOError:
        return False
    if ensembler_stage in project_metadata.data:
        return True
    else:
        return False


def check_ensembler_modeling_stage_complete(ensembler_stage, targetid):
    if not check_ensembler_modeling_stage_metadata_exists(ensembler_stage, targetid):
        return False
    if not check_ensembler_modeling_stage_first_model_file_exists(ensembler_stage, targetid):
        return False
    return True


def select_templates_by_seqid_cutoff(targetid, seqid_cutoff=None):
    """
    :param seqid_cutoff:
    :return:
    """
    seqid_filepath = os.path.join(ensembler.core.default_project_dirnames.models, targetid, 'sequence-identities.txt')
    with open(seqid_filepath) as seqid_file:
        seqid_lines_split = [line.split() for line in seqid_file.read().splitlines()]

    templateids = np.array([i[0] for i in seqid_lines_split])
    seqids = np.array([float(i[1]) for i in seqid_lines_split])

    # must coerce to string due to yaml.dump type requirements
    selected_templateids = [str(x) for x in templateids[seqids > seqid_cutoff]]

    return selected_templateids