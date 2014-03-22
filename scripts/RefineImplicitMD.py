#!/usr/bin/env python
#
# Refine models with implicit-solvent MD simulations
#
# Daniel L. Parton <partond@mskcc.org> - 21 Mar 2014
#

import MSMSeeder
import MSMSeeder.refinement
import mpi4py.MPI
comm = mpi4py.MPI.COMM_WORLD 
rank = comm.rank
size = comm.size

# ========
# Parse command-line arguments
# ========

import argparse
argparser = argparse.ArgumentParser(description='Models a set of target sequences onto a set of template structures using Modeller.', formatter_class=argparse.RawTextHelpFormatter)

argparser.add_argument('--OpenMMPlatform', choices=['CUDA', 'OpenCL', 'CPU', 'Reference'], default='CUDA', help='(Default: CUDA) Choose the OpenMM Platform to use.')
argparser.add_argument('-gpupn', type=int, default=1, help='(Default: 1) If using GPUs, select how many are available per node.')
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
# Run simulations
# ========

MSMSeeder.refinement.refine_implicitMD(openmm_platform=args.OpenMMPlatform, gpupn=args.gpupn, process_only_these_targets=args.ProcessOnlyTheseTargets, process_only_these_templates=args.ProcessOnlyTheseTemplates)

