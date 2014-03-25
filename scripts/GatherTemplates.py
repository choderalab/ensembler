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
database, specified via either the project metadata file or the --db_path
argument.

"UniProt": Select templates from UniProt with a user-defined query
string, plus an optional subquery regex string.'''

helpstring_uniprot_query = (
r'''Specify a UniProt search string, using the same syntax as on the UniProt
site (note: not all syntax may be supported, but most basic searches will
work). *All* domains contained within the returned UniProt entries will be
selected as templates, unless the --uniprot_domain_regex option is used to
select a subset. The script will print some information on the set of unique
domain names returned by the initial UniProt search, which can help with
constructing a suitable string for --uniprot_domain_regex.

Example: 'domain:"Protein kinase" AND reviewed:yes' - this will return reviewed
UniProt entries containing "Protein kinase" domain annotations, but note that
all domains contained with those entries (including domains which are not
"Protein kinase") will be selected as templates, unless the
--uniprot_domain_refex flag is also set.
''')

helpstring_uniprot_domain_regex = (
r'''Optional regular expression for subselecting domains from within UniProt
entries (case-sensitive). If not provided, all domains contained within
returned UniProt entries will be selected as templates (this will often not be
the desired behavior).

Example: '^Protein kinase(?!; truncated)(?!; inactive)' - matches "Protein
kinase" as well as "Protein kinase; 1" and "Protein kinase; 2"
''')

argparser.add_argument('--gather_from', type=str, help=helpstring_gatherfrom, choices=['TargetExplorerDB', 'UniProt'])
argparser.add_argument('--db_path', type=str, help='TargetExplorerDB database path', default=None)
argparser.add_argument('--uniprot_query', type=str, help=helpstring_uniprot_query, default=None)
argparser.add_argument('--uniprot_domain_regex', type=str, help=helpstring_uniprot_domain_regex, default=None)
argparser.add_argument('--structure_paths', type=list, help='Optional list of local directories within which to search for PDB and SIFTS files.', default=None)
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
template_selection_method = args.gather_from

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
    DB_path = args.db_path

    # Otherwise check project metadata file
    if DB_path == None:
        DB_path = project_metadata.get(('template-selection', 'TargetExplorer-database-path'))

    if DB_path == None:
        raise Exception, 'Database path not found in command-line args or in project metadata file. Cannot continue'

elif template_selection_method == 'UniProt':

    # Command-line args take priority
    UniProt_query_string = args.uniprot_query
    UniProt_domain_regex = args.uniprot_domain_regex
    structure_paths = args.structure_paths

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

