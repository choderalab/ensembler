#!/usr/bin/env python
#
# Refine models with implicit-solvent MD simulations
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 24 Mar 2014
#

import argparse
import msmseeder
import msmseeder.packaging

def main():
    # ========
    # Parse command-line arguments
    # ========

    argparser = argparse.ArgumentParser(description='Packages models for transfer or for set-up as a Folding@Home project. Packaged models are stored in the packaged-models/ directory.', formatter_class=argparse.RawTextHelpFormatter)

    package_for_helpstring = r'''
    Specify which packaging method to use.

    transfer: compress results into a single .tgz file
    FAH: set-up the input files and directory structure necessary to start a Folding@Home project.
    '''.strip()

    argparser.add_argument('--package_for', choices=['transfer', 'FAH'], help=package_for_helpstring)
    argparser.add_argument('--nFAHclones', type=int, default=10, help='(Default: 10) If packaging for Folding@Home, select the number of clones to use for each model.')
    argparser.add_argument('--archiveFAHproject', type=bool, default=True, help='(Default: True) If packaging for Folding@Home, choose whether to compress the results into a .tgz file.')
    argparser.add_argument('--targets', nargs='+', help='(Default: all targets) Optionally define a subset of targets to work on by providing one or more target IDs separated by spaces (e.g. "ABL1_HUMAN_D0")')
    args = argparser.parse_args()

    msmseeder.core.check_project_toplevel_dir()

    # ========
    # Run simulations
    # ========

    if args.package_for == 'transfer':
        msmseeder.packaging.package_for_transfer(process_only_these_targets=args.targets)

    elif args.package_for == 'FAH':
        msmseeder.packaging.package_for_fah(process_only_these_targets=args.targets, nclones=args.nFAHclones, archive=args.archiveFAHproject)

if __name__ == '__main__':
    main()
