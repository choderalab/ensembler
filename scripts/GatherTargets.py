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
database, specified via the --db_path
argument.

"UniProt": Gather target data from UniProt with a user-defined query
string.'''
argparser.add_argument('--gather_from', type=str, help=helpstring_gatherfrom, choices=['TargetExplorerDB', 'UniProt'], default='TargetExplorerDB')
argparser.add_argument('--db_path', type=str, help='TargetExplorerDB database path. Will be converted to an absolute path.', required=True)
argparser.add_argument('--species', type=str, help='Optional NCBI taxonomy ID specifier e.g. 9606 (human).')
args = argparser.parse_args()

msmseeder.core.check_project_toplevel_dir()

# ========
# Get the target selection method
# ========

target_selection_method = args.gather_from

# ========
# Get method-specific parameters
# ========

if target_selection_method == 'TargetExplorerDB':
    DB_path = args.db_path

elif target_selection_method == 'UniProt':
    pass


# ========
# Run the selected gather targets method
# ========

if target_selection_method == 'TargetExplorerDB':
    msmseeder.initproject.gather_targets_from_TargetExplorerDB(DB_path, species=args.species)

elif target_selection_method == 'UniProt':
    msmseeder.initproject.gather_targets_from_UniProt()

