# ========
# Global package variables
# ========

import os
import msmseeder

src_toplevel_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

datestamp_format_string = '%Y-%m-%d %H:%M:%S UTC'

project_metadata_filename = 'project-data.yaml'
manual_overrides_filename = 'manual-overrides.yaml'

template_acceptable_ratio_observed_residues = 0.7

# listed in order
msmseeder_stages = [
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

project_dirnames = [
    'targets', 'structures', 'templates', 'models', 'packaged-models',
    os.path.join('structures', 'pdb'),
    os.path.join('structures', 'sifts'),
    os.path.join('templates', 'structures')
]

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


class ProjectMetadata:
    def __init__(self, data):
        self.data = data

    def write(self, ofilepath):
        with open(ofilepath, 'w') as ofile:
            for stage in msmseeder_stages:
                if stage in self.data.keys():
                    subdict = {stage: self.data[stage]}
                    yaml.dump(subdict, ofile, default_flow_style=False)


def write_metadata(new_metadata_dict, msmseeder_stage):
    if msmseeder_stage == 'init':
        metadata_dict = {}
    else:
        prev_msmseeder_stage = msmseeder_stages[msmseeder_stages.index(msmseeder_stage) - 1]
        prev_metadata_filepath = metadata_file_mapper(prev_msmseeder_stage)
        with open(prev_metadata_filepath) as prev_metadata_file:
            metadata_dict = yaml.load(prev_metadata_file)

    metadata_dict.update(new_metadata_dict)
    metadata = ProjectMetadata(metadata_dict)
    metadata.write(metadata_file_mapper(msmseeder_stage))


def metadata_file_mapper(msmseeder_stage, target_id=None):
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
    target_fasta_string = '>%s\n%s\n' % (seqid, seq)
    return target_fasta_string