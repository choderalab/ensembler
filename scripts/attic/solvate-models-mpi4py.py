# Subject models to explicit solvent simulation.
#
# John D. Chodera <choderaj@mskcc.org> - 17 Feb 2013
#
# PREREQUISITES
#
# * OpenMM
# http://simtk.org/home/openmm
#
# TODO
# * use ff99sb-ildn-star
# * trim all systems to have the same number of waters?

# PARAMETERS

import sys
from ast import literal_eval

# Process only these targets, if specified.
# e.g. -targets '["SRC_HUMAN_PK0_P12931", "ABL1_HUMAN_PK0_P00519"]'
try:
    process_only_these_targets = literal_eval( sys.argv[ sys.argv.index('-targets') + 1 ] )
except ValueError:
    process_only_these_targets = False

if process_only_these_targets:
    print 'Processing only these targets:'
    print process_only_these_targets

# OpenMM parameters

import simtk.openmm as mm
import simtk.unit as units
import simtk.openmm.app as app

forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml'] # list of forcefields to use in parameterization
nparticles_per_water = 3 # number of particles per water molecule

#box_width = 90.0 * units.angstroms
#boxsize = box_width * mm.Vec3(1,1,1)
padding = 10.0 * units.angstroms

verbose = True
write_trajectory = False

#
# GET ABSOLUTE PATHS
#

import os.path

# Input files.
targets_directory = os.path.abspath("targets") # target sequences for modeling
templates_directory = os.path.abspath("templates") # template structures for use in modeling
models_directory = os.path.abspath("models")

#
# READ TEMPLATE AND TARGET INDICES
#

targets_index_filename = os.path.join(targets_directory, 'targets.txt')
infile = open(targets_index_filename, 'r')
targets = [ line.strip() for line in infile ]
infile.close()
print '%d target sequences' % len(targets)
#print "targets:"
#print targets

templates_index_filename = os.path.join(templates_directory, 'templates.txt')
infile = open(templates_index_filename, 'r')
templates = [ line.strip() for line in infile ]
infile.close()
print '%d template structures' % len(templates)
#print "templates:"
#print templates

#
# LOAD FORCEFIELD
#

forcefield = app.ForceField(*forcefields_to_use)

#
# MPI
#

import mpi4py
from mpi4py import MPI

# get communicator
comm = MPI.COMM_WORLD 

# get rank and size
rank = comm.rank
size = comm.size

#
# SOLVATE MODELS
#

original_directory = os.getcwd()

for target in targets:
    
    # Process only specified targets if directed.
    if process_only_these_targets and (target not in process_only_these_targets): continue

    target_directory = os.path.join(models_directory, target)
    if not os.path.exists(target_directory): continue

    # Start a 'reject file'.
    reject_filename = os.path.join(target_directory, 'reject-explicit.txt')
    reject_file = open(reject_filename, 'w')

    # Process all templates.
    for template_index in range(rank, len(templates), size):
        template = templates[template_index]

        print "-------------------------------------------------------------------------"
        print "Solvating %s => %s in explicit solvent" % (target, template)
        print "-------------------------------------------------------------------------"
        
        model_directory = os.path.join(models_directory, target, template)
        if not os.path.exists(model_directory): continue

        os.chdir(model_directory)

        model_filename = os.path.join(model_directory, 'implicit-refined.pdb')
        if not os.path.exists(model_filename): continue


	# Pass if this simulation has already been run.
        nwaters_filename = os.path.join(model_directory, 'nwaters.txt')
	if os.path.exists(nwaters_filename): continue

        try:
            if verbose: print "Reading model..."
            pdb = app.PDBFile(model_filename)

            # Count initial atoms.
            natoms_initial = len(pdb.positions)

            if verbose: print "Solvating model..."
            modeller = app.Modeller(pdb.topology, pdb.positions)
            #modeller.addSolvent(forcefield, model='tip3p', boxSize=boxsize)
            modeller.addSolvent(forcefield, model='tip3p', padding=padding)
            topology = modeller.getTopology()
            positions = modeller.getPositions()

            # Count final atoms.
            natoms_final = len(positions)
            nwaters = (natoms_final - natoms_initial) / nparticles_per_water
            if verbose: print "Solvated model contains %d waters" % nwaters

            # Record waters.
            outfile = open(nwaters_filename, 'w')
            outfile.write('%d\n' % nwaters)
            outfile.close()

            os.chdir(original_directory)    

        except Exception as e:
            # Add to rejection file.
            reject_file.write('%s : %s\n' % (template, str(e)))
            reject_file.flush()
            

