#!/usr/bin/env python
#
# Refine models with implicit-solvent MD simulations
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

    argparser = argparse.ArgumentParser(description='Conducts implicit-solvent MD refinement on a set of models.', formatter_class=argparse.RawTextHelpFormatter)

    argparser.add_argument('--openmm_platform', choices=['CUDA', 'OpenCL', 'CPU', 'Reference'], default='CUDA', help='(Default: CUDA) Choose the OpenMM Platform to use.')
    argparser.add_argument('-gpupn', type=int, default=1, help='(Default: 1) If using GPUs, select how many are available per node.')
    argparser.add_argument('--targets', nargs='+', help='(Default: all targets) Optionally define a subset of targets to work on by providing one or more target IDs separated by spaces (e.g. "ABL1_HUMAN_D0")')
    argparser.add_argument('--templates', nargs='+', help='(Default: all templates) Optionally define a subset of templates to work on by providing one or more template IDs separated by spaces (e.g. "ABL1_HUMAN_D0_1OPL_A")')
    argparser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    args = argparser.parse_args()

    msmseeder.core.check_project_toplevel_dir()

    # ========
    # Run simulations
    # ========

    msmseeder.refinement.refine_implicitMD(openmm_platform=args.openmm_platform, gpupn=args.gpupn, process_only_these_targets=args.targets, process_only_these_templates=args.templates, verbose=args.verbose)

if __name__ == '__main__':
    main()
