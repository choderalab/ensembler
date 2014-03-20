#!/usr/bin/env python
#
# Models a set of targets sequences onto a set of template structures using Modeller.
#
# Daniel L. Parton <partond@mskcc.org> - 11 Mar 2014
#

import MSMSeeder
import MSMSeeder.modelling

# ========
# Parse command-line arguments
# ========

import argparse
argparser = argparse.ArgumentParser(description='Models a set of target sequences onto a set of template structures using Modeller.', formatter_class=argparse.RawTextHelpFormatter)

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
# Run the selected gather templates method
# ========

MSMSeeder.modelling.build_models(process_only_these_targets=args.ProcessOnlyTheseTargets, process_only_these_templates=args.ProcessOnlyTheseTemplates)
