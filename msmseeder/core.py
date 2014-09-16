# ========
# Global package variables
# ========

import os
src_toplevel_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

datestamp_format_string = '%Y-%m-%d %H:%M:%S UTC'

project_metadata_filename = 'project-data.yaml'
manual_overrides_filename = 'manual-overrides.yaml'

template_acceptable_ratio_observed_residues = 0.7

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
    def __init__(self, log_filepath, additional_log_data={}):
        import socket
        self.log_filepath = log_filepath
        log_data = {
            'datestamp': get_utcnow_formatted(),
            'hostname': socket.gethostname(),
        }

        log_data.update(additional_log_data)

        with open(self.log_filepath, 'w') as log_file:
            yaml.dump(log_data, log_file, default_flow_style=False)

    def log(self, new_log_data):
        with open(self.log_filepath, 'r') as log_file:
            log_data = yaml.load(log_file)

        log_data.update(new_log_data)

        with open(self.log_filepath, 'w') as log_file:
            yaml.dump(log_data, log_file, default_flow_style=False)

class ProjectMetadata:
    def __init__(self, data):
        # Listed in desired document order
        self.project_metadata_categories = ['init', 'gather_targets', 'gather_templates', 'build_models', 'sort_by_sequence_identity', 'cluster_models', 'refine_implicit_md', 'solvate_models', 'determine_nwaters', 'refine_explicit_md', 'package_for_fah']
        self.data = data

    def write(self, ofilepath):
        with open(ofilepath, 'w') as ofile:
            for category in self.project_metadata_categories:
                if category in self.data.keys():
                    subdict = {category: self.data[category]}
                    yaml.dump(subdict, ofile, default_flow_style=False)

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
