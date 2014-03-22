#!/usr/bin/env python
#
# Solvate models which have been through implicit-solvent MD equilibration
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 21 Mar 2014
#

import MSMSeeder
import MSMSeeder.refinement

# ========
# Parse command-line arguments
# ========

import argparse
argparser = argparse.ArgumentParser(description='Solvates each model individually, then determines the distribution of the number of waters in each model and selects the value at the 68th percentile.', formatter_class=argparse.RawTextHelpFormatter)

argparser.add_argument('--ProcessOnlyTheseTargets', nargs='+', help='Supply one or more target IDs separated by spaces (e.g. "ABL1_HUMAN_D0")')
argparser.add_argument('--ProcessOnlyTheseTemplates', nargs='+', help='Supply one or more template IDs separated by spaces (e.g. "ABL1_HUMAN_D0_1OPL_A")')
args = argparser.parse_args()

MSMSeeder.core.check_project_toplevel_dir()

# ========
# Parse project metadata
# ========

project_metadata = MSMSeeder.core.ProjectMetadata()
project_metadata.load(MSMSeeder.core.project_metadata_filename)

# ========
# Solvate each model individually
# ========

MSMSeeder.refinement.solvate_models(process_only_these_targets=args.ProcessOnlyTheseTargets, process_only_these_templates=args.ProcessOnlyTheseTemplates)

# ========
# Determine distribution of nwaters in each model, and select the value at the 68th percentile
# ========

MSMSeeder.refinement.determine_nwaters(process_only_these_targets=args.ProcessOnlyTheseTargets, process_only_these_templates=args.ProcessOnlyTheseTemplates)

