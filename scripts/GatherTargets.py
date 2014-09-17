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
database, specified via the --dbapi_uri
argument.

"UniProt": Gather target data from UniProt with a user-defined query
string.'''
argparser.add_argument('--gather_from', type=str, help=helpstring_gatherfrom, choices=['TargetExplorerDB', 'UniProt'], default='TargetExplorerDB')
argparser.add_argument('--dbapi_uri', type=str, help='TargetExplorer database API URI, e.g. "http://plfah2.mskcc.org/kinomeDBAPI"', required=True)
argparser.add_argument('--query', type=str, help='Query string for TargetExplorer database API, e.g. \'species="Human"\'')
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
    msmseeder.initproject.gather_targets_from_uniprot()
