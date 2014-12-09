#!/usr/bin/env python
#
# Solvate models which have been through implicit-solvent MD equilibration
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 21 Mar 2014
#

import argparse
import msmseeder
import msmseeder.refinement

def main():
    # ========
    # Parse command-line arguments
    # ========

    argparser = argparse.ArgumentParser(description='Solvates each model individually, then determines the distribution of the number of waters in each model and selects the value at the 68th percentile.', formatter_class=argparse.RawTextHelpFormatter)

    argparser.add_argument('--targets', nargs='+', help='(Default: all targets) Optionally define a subset of targets to work on by providing one or more target IDs separated by spaces (e.g. "ABL1_HUMAN_D0")')
    argparser.add_argument('--templates', nargs='+', help='(Default: all templates) Optionally define a subset of templates to work on by providing one or more template IDs separated by spaces (e.g. "ABL1_HUMAN_D0_1OPL_A")')
    args = argparser.parse_args()

    msmseeder.core.check_project_toplevel_dir()

    # ========
    # Solvate each model individually
    # ========

    msmseeder.refinement.solvate_models(process_only_these_targets=args.targets, process_only_these_templates=args.templates)

    # ========
    # Determine distribution of nwaters in each model, and select the value at the 68th percentile
    # ========

    msmseeder.refinement.determine_nwaters(process_only_these_targets=args.targets, process_only_these_templates=args.templates)

if __name__ == '__main__':
    main()
