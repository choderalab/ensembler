#!/usr/bin/env python
#
# Models a set of targets sequences onto a set of template structures using Modeller.
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 11 Mar 2014
#

import argparse
import msmseeder
import msmseeder.modelling

def main():
    argparser = argparse.ArgumentParser(description='Sorts models by sequence identity, then performs clustering to filter out non-unique models.', formatter_class=argparse.RawTextHelpFormatter)
    argparser.add_argument('--targets', nargs='+', help='(Default: all targets) Optionally define a subset of targets to work on by providing one or more target IDs separated by spaces (e.g. "ABL1_HUMAN_D0")')
    argparser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    args = argparser.parse_args()

    msmseeder.core.check_project_toplevel_dir()

    # if args.verbose:
    #     loglevel = 'debug'
    # else:
    #     loglevel = 'info'

    msmseeder.modelling.sort_by_sequence_identity(process_only_these_targets=args.targets, loglevel=None)
    msmseeder.modelling.cluster_models(process_only_these_targets=args.targets, verbose=args.verbose)

if __name__ == '__main__':
    main()
