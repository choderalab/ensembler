#!/usr/bin/env python
#
# Gathers protein target data - IDs and sequences.
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 11 Mar 2014
#

import msmseeder
import msmseeder.initproject

# ========
# Parse command-line arguments
# ========

import argparse
argparser = argparse.ArgumentParser(description='Gather target protein data - IDs and sequences.', formatter_class=argparse.RawTextHelpFormatter)

helpstring_gatherfrom = r'''Choose a method for gathering target data.

"TargetExplorerDB": Gather target data from an existing TargetExplorer
database, specified via the --dbapi_uri argument.

"UniProt": Select targets from UniProt with a user-defined query
string, plus an optional subquery regex string.'''

helpstring_uniprot_query = (
r'''Specify a UniProt search string, using the same syntax as on the UniProt
site (note: not all syntax may be supported, but most basic searches will
work). *All* domains contained within the returned UniProt entries will be
selected as targets, unless the --uniprot_domain_regex option is used to
select a subset. The script will print some information on the set of unique
domain names returned by the initial UniProt search, which can help with
constructing a suitable string for --uniprot_domain_regex.

Example: 'domain:"Protein kinase" AND taxonomy:9606 AND reviewed:yes' - this
will return reviewed UniProt entries for human (taxonomy ID: 9606) proteins
containing "Protein kinase" domain annotations. Note that all domains contained
with those entries (including domains which are not "Protein kinase") will be
selected as targets, unless the --uniprot_domain_regex flag is also set.
''')

helpstring_uniprot_domain_regex = (
r'''Optional regular expression for subselecting domains from within UniProt
entries (case-sensitive). If not provided, all domains contained within
returned UniProt entries will be selected as targets (this will often not be
the desired behavior).

Example: '^Protein kinase(?!; truncated)(?!; inactive)' - matches "Protein
kinase" as well as "Protein kinase; 1" and "Protein kinase; 2"
''')

argparser.add_argument('--gather_from', type=str, help=helpstring_gatherfrom, choices=['TargetExplorerDB', 'UniProt'], default='TargetExplorerDB')
argparser.add_argument('--dbapi_uri', type=str, help='TargetExplorer database API URI, e.g. "http://plfah2.mskcc.org/kinomeDBAPI"')
argparser.add_argument('--query', type=str, help='Query string for TargetExplorer database API, e.g. \'species="Human"\'')
argparser.add_argument('--uniprot_query', type=str, help=helpstring_uniprot_query)
argparser.add_argument('--uniprot_domain_regex', type=str, help=helpstring_uniprot_domain_regex)
args = argparser.parse_args()

msmseeder.core.check_project_toplevel_dir()

# ========
# Get the target selection method
# ========

target_selection_method = args.gather_from

# ========
# Run the selected gather targets method
# ========

if target_selection_method == 'TargetExplorerDB':
    msmseeder.initproject.gather_targets_from_targetexplorerdb(args.dbapi_uri, search_string=args.query)

elif target_selection_method == 'UniProt':
    msmseeder.initproject.gather_targets_from_uniprot(args.uniprot_query, args.uniprot_domain_regex)
