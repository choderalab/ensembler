#!/usr/bin/env python
#
# Gathers protein template data - IDs, sequences and structures.
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 11 Mar 2014
#

import os
import msmseeder
import msmseeder.initproject

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

argparser.add_argument('--gather_from', type=str, help=helpstring_gatherfrom, choices=['TargetExplorerDB', 'UniProt'], default='UniProt')
argparser.add_argument('--db_path', type=str, help='TargetExplorerDB database path', default=None)
argparser.add_argument('--uniprot_query', type=str, help=helpstring_uniprot_query, default=None)
argparser.add_argument('--uniprot_domain_regex', type=str, help=helpstring_uniprot_domain_regex, required=True)
argparser.add_argument('--structure_paths', help='(Optional) Local directories within which to search for PDB and SIFTS files (space-separated).', nargs='*')
args = argparser.parse_args()

msmseeder.core.check_project_toplevel_dir()

# ========
# Get the template selection method
# ========

template_selection_method = args.gather_from

# ========
# Get method-specific parameters
# ========

if template_selection_method == 'TargetExplorerDB':

    DB_path = args.db_path

elif template_selection_method == 'UniProt':

    UniProt_query_string = args.uniprot_query
    UniProt_domain_regex = args.uniprot_domain_regex
    if args.structure_paths != None:
        structure_paths = [os.path.abspath(path) for path in args.structure_paths]
    else:
        structure_paths = []

# ========
# Run the selected gather templates method
# ========

if template_selection_method == 'TargetExplorerDB':
    msmseeder.initproject.gather_templates_from_TargetExplorerDB(DB_path)

if template_selection_method == 'UniProt':
    msmseeder.initproject.gather_templates_from_UniProt(UniProt_query_string, UniProt_domain_regex, structure_paths)

