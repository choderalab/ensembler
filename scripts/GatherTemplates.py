#!/usr/bin/env python
#
# Gathers protein template data - IDs, sequences and structures.
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 11 Mar 2014
#

import MSMSeeder
import MSMSeeder.initproject

# ========
# Parse command-line arguments
# ========

import argparse
argparser = argparse.ArgumentParser(description='Gather template protein data - IDs, sequences and structures.', formatter_class=argparse.RawTextHelpFormatter)

helpstring_gatherfrom = r'''Choose a method for selecting the templates.

"TargetExplorerDB": Select templates from an existing TargetExplorerDB
database, specified via either the project metadata file or the --DBpath
argument.

"UniProt": Select templates from UniProt with a user-defined query
string.'''
argparser.add_argument('--GatherFrom', type=str, help=helpstring_gatherfrom, choices=['TargetExplorerDB', 'UniProt'])
argparser.add_argument('--DBpath', type=str, help='TargetExplorerDB database path', default=None)
argparser.add_argument('--UniProtQuery', type=str, help='UniProt query string', default=None)
argparser.add_argument('--UniProtDomainRegex', type=str, help='Optional regular expression for selecting domains from within UniProt entries (case-sensitive)', default=None)
argparser.add_argument('--StructurePaths', type=list, help='Optional list of local directories within which to search for PDB and SIFTS files.', default=None)
args = argparser.parse_args()

MSMSeeder.core.check_project_toplevel_dir()

# ========
# Parse project metadata
# ========

project_metadata = MSMSeeder.core.ProjectMetadata()
project_metadata.load(MSMSeeder.core.project_metadata_filename)

# ========
# Get the template selection method
# ========

# Command-line args take priority
template_selection_method = args.GatherFrom

# Otherwise check project metadata file
if template_selection_method == None:
    template_selection_method = project_metadata.get(('template-selection', 'template-selection-method'))

if template_selection_method == None:
    raise Exception, 'Template selection method not found in command-line args or in project metadata file. Cannot continue'

# ========
# Get method-specific parameters
# ========

if template_selection_method == 'TargetExplorerDB':

    # Command-line args take priority
    DB_path = args.DBpath

    # Otherwise check project metadata file
    if DB_path == None:
        DB_path = project_metadata.get(('template-selection', 'TargetExplorer-database-path'))

    if DB_path == None:
        raise Exception, 'Database path not found in command-line args or in project metadata file. Cannot continue'

elif template_selection_method == 'UniProt':

    # Command-line args take priority
    UniProt_query_string = args.UniProtQuery
    UniProt_domain_regex = args.UniProtDomainRegex
    structure_paths = args.StructurePaths

    # Otherwise check project metadata file
    if UniProt_query_string == None:
        UniProt_query_string = project_metadata.get(('template-selection', 'UniProt-query-string'))
    if UniProt_query_string == None:
        raise Exception, 'UniProt query string not found in command-line args or in project metadata file. Cannot continue'

    if UniProt_domain_regex == None:
        UniProt_domain_regex = project_metadata.get(('template-selection', 'UniProt-domain-regex'))

    if structure_paths == None:
        structure_paths = project_metadata.get(('template-selection', 'structure-paths'))
    if structure_paths == None:
        structure_paths = []

# ========
# Run the selected gather templates method
# ========

if template_selection_method == 'TargetExplorerDB':
    MSMSeeder.initproject.gather_templates_from_TargetExplorerDB(DB_path)

if template_selection_method == 'UniProt':
    MSMSeeder.initproject.gather_templates_from_UniProt(UniProt_query_string, UniProt_domain_regex, structure_paths)

