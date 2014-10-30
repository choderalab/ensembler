import os
import logging
import sys
import re
from collections import namedtuple

# ========
# Global package variables
# ========
import Bio

src_toplevel_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

datestamp_format_string = '%Y-%m-%d %H:%M:%S UTC'

project_metadata_filename = 'project-data.yaml'
manual_overrides_filename = 'manual-overrides.yaml'

template_acceptable_ratio_observed_residues = 0.7

# listed in order
project_stages = [
    'init',
    'gather_targets',
    'gather_templates',
    'build_models',
    'sort_by_sequence_identity',
    'cluster_models',
    'refine_implicit_md',
    'solvate_models',
    'determine_nwaters',
    'refine_explicit_md',
    'package_for_fah',
]

ProjectDirNames = namedtuple(
    'ProjectDirNames',
    ['targets', 'templates', 'structures', 'models', 'packaged_models', 'structures_pdb', 'structures_sifts', 'templates_structures']
)

default_project_dirnames = ProjectDirNames(
    targets='targets',
    templates='templates',
    structures='structures',
    models='models',
    packaged_models='packaged_models',
    structures_pdb=os.path.join('structures', 'pdb'),
    structures_sifts=os.path.join('structures', 'sifts'),
    templates_structures=os.path.join('templates', 'structures'),
)

logger = logging.getLogger('info')
default_loglevel = 'info'
loglevel_obj = getattr(logging, default_loglevel.upper())
logger.setLevel(loglevel_obj)
logger.addHandler(logging.StreamHandler(stream=sys.stdout))

metadata_filename_regex = re.compile('(meta)([0-9]+)\.yaml')

# ========
# YAML
# ========

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

def check_project_toplevel_dir():
    import os
    for dirpath in ['structures', 'templates', 'targets', 'models']:
        if not os.path.exists(dirpath):
            raise Exception, 'Current directory not recognized as the top-level directory of a project.'

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
            yaml.dump(self.log_data, log_file, default_flow_style=False)


class ManualOverrides:
    def __init__(self):
        import yaml
        if os.path.exists(manual_overrides_filename):
            with open(manual_overrides_filename, 'r') as manual_overrides_file:
                manual_overrides_yaml = yaml.load(manual_overrides_file)
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


def metadata_dir_mapper(msmseeder_stage, target_id=None):
    metadata_dir_dict = {
        'init': '.',
        'gather_targets': 'targets',
        'gather_templates': 'templates',
    }
    if msmseeder_stage in metadata_dir_dict:
        return metadata_dir_dict[msmseeder_stage]
    elif msmseeder_stage in ['build_models', 'sort_by_sequence_identity', 'cluster_models', 'refine_implicit_md', 'solvate_models', 'determine_nwaters', 'refine_explicit_md']:
        return os.path.join('models', target_id)


class ProjectMetadata:
    """
    Examples
    --------
    >>> init_project_metadata = msmseeder.core.ProjectMetadata(project_stage='init')
    >>> test_data = {'test_field': 'test_value'}
    >>> init_project_metadata.add_data(test_data)
    >>> init_project_metadata.write()
    >>> gather_targets_project_metadata = msmseeder.core.ProjectMetadata(project_stage='gather_targets')
    >>> gather_targets_project_metadata.add_data(test_data)
    >>> gather_targets_project_metadata.write()
    """
    def __init__(self, project_stage='init'):
        self.data = {}
        self.project_stage = project_stage
        if project_stage != project_stages[0]:
            self.add_all_prev_metadata(project_stage)

    def add_all_prev_metadata(self, project_stage):
        current_project_stage_index = project_stages.index(project_stage)
        for prev_project_stage in project_stages[: current_project_stage_index]:
            self.add_prev_metadata(prev_project_stage)

    def add_prev_metadata(self, project_stage):
        latest_metadata_filepath = self.determine_latest_metadata_file(project_stage)
        with open(latest_metadata_filepath) as latest_metadata_file:
            prev_metadata = yaml.load(latest_metadata_file)
        self.add_data(prev_metadata[project_stage], project_stage=project_stage)

    def determine_latest_metadata_file(self, project_stage):
        metadata_dir = metadata_dir_mapper(project_stage)
        latest_metadata_file_index = self.determine_latest_metadata_file_index(project_stage)
        latest_metadata_filepath = self.gen_metadata_filepath_from_dir_and_index(metadata_dir, latest_metadata_file_index)
        return latest_metadata_filepath

    def determine_latest_metadata_file_index(self, project_stage):
        """
        Returns -1 if no metadata files found
        :param project_stage: str
        :return: int
        """
        metadata_dir = metadata_dir_mapper(project_stage)
        dir_contents = os.listdir(metadata_dir)
        metadata_file_indices = []
        for filename in dir_contents:
            match = re.match(metadata_filename_regex, filename)
            if match:
                metadata_file_indices.append(int(match.groups()[1]))
        if len(metadata_file_indices) > 0:
            return max(metadata_file_indices)
        else:
            return -1

    def gen_metadata_filepath_from_dir_and_index(self, dirpath, index):
        metadata_filepath = os.path.join(dirpath, 'meta%d.yaml' % index)
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
        >>> project_metadata = msmseeder.core.ProjectMetadata(project_stage='init')
        >>> metadata = {'datestamp': msmseeder.core.get_utcnow_formatted()}
        >>> project_metadata.add_data(metadata)
        """
        if project_stage is None:
            project_stage = self.project_stage
        self.data.update({
            project_stage: data
        })

    def write(self):
        metadata_dir = metadata_dir_mapper(self.project_stage)
        latest_metadata_file_index = self.determine_latest_metadata_file_index(self.project_stage)
        metadata_filepath = self.gen_metadata_filepath_from_dir_and_index(metadata_dir, latest_metadata_file_index+1)
        with open(metadata_filepath, 'w') as ofile:
            for stage in project_stages:
                if stage in self.data.keys():
                    subdict = {stage: self.data[stage]}
                    yaml.dump(subdict, ofile, default_flow_style=False)


class DeprecatedProjectMetadata:
    # TODO deprecate
    def __init__(self, data):
        self.data = data

    def write(self, ofilepath):
        with open(ofilepath, 'w') as ofile:
            for stage in project_stages:
                if stage in self.data.keys():
                    subdict = {stage: self.data[stage]}
                    yaml.dump(subdict, ofile, default_flow_style=False)


def write_metadata(new_metadata_dict, msmseeder_stage, target_id=None):
    # TODO deprecate
    if msmseeder_stage == 'init':
        metadata_dict = {}
    else:
        prev_msmseeder_stage = project_stages[project_stages.index(msmseeder_stage) - 1]
        prev_metadata_filepath = metadata_file_mapper(prev_msmseeder_stage, target_id=target_id)
        with open(prev_metadata_filepath) as prev_metadata_file:
            metadata_dict = yaml.load(prev_metadata_file)

    metadata_dict.update(new_metadata_dict)
    metadata = ProjectMetadata(metadata_dict)
    metadata.write(metadata_file_mapper(msmseeder_stage, target_id=target_id))


def metadata_file_mapper(msmseeder_stage, target_id=None):
    # TODO deprecate
    metadata_file_dict = {
        'init': 'meta.yaml',
        'gather_targets': os.path.join('targets', 'meta.yaml'),
        'gather_templates': os.path.join('templates', 'meta.yaml'),
    }
    if msmseeder_stage in metadata_file_dict:
        return metadata_file_dict[msmseeder_stage]
    elif msmseeder_stage in ['build_models', 'sort_by_sequence_identity', 'cluster_models', 'refine_implicit_md', 'solvate_models', 'determine_nwaters', 'refine_explicit_md']:
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


def construct_fasta_str(seqid, seq):
    target_fasta_string = '>%s\n%s\n' % (seqid, seqwrap(seq).strip())
    return target_fasta_string


def get_targets():
    targets_dir = os.path.abspath('targets')
    targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
    targets = list(Bio.SeqIO.parse(targets_fasta_filename, 'fasta'))
    return targets


def get_templates():
    templates_dir = os.path.abspath('templates')
    templates_fasta_filename = os.path.join(templates_dir, 'templates.fa')
    templates = list(Bio.SeqIO.parse(templates_fasta_filename, 'fasta'))
    return templates


def get_targets_and_templates():
    targets = get_targets()
    templates = get_templates()
    return targets, templates