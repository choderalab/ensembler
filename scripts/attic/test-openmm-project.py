# Test to make sure CLONEs were packaged correctly using short MD simulations.
#
# John D. Chodera <choderaj@mskcc.org> - 27 Mar 2013
#
# PREREQUISITES
#
# * OpenMM
# http://simtk.org/home/openmm

# PARAMETERS

import sys
from ast import literal_eval

# XXX NOTE: process_only_these_targets is not used in this script. Instead modify variable: project_to_test
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

nsteps = 10 # number of steps to run

platform_name = "CUDA"
try:
        ngpus = int( sys.argv[ sys.argv.index('-gpupn') + 1 ] )
except ValueError:
        ngpus = 1 # number of GPUs per node

verbose = True

#
# GET ABSOLUTE PATHS
#

import os.path


projects_directory = os.path.abspath("projects")
project_to_test = "ABL1_HUMAN_PK0_P00519"

project_directory = os.path.join(projects_directory, project_to_test)

#
# SUBROUTINES
#

def read_file(filename):
    infile = open(filename, 'r')
    contents = infile.read()
    infile.close()
    return contents

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
gpuid = (rank % ngpus)

# Choose platform.
platform = mm.Platform.getPlatformByName(platform_name)

# Set GPU id.
platform.setPropertyDefaultValue('CudaDeviceIndex', '%d' % gpuid)
platform.setPropertyDefaultValue('OpenCLDeviceIndex', '%d' % gpuid)

# Report
import socket
print "host: %s, MPI rank %d / %d, gpuid %d" % (socket.gethostname(), rank, size, gpuid)

#
# SIMULATE MODELS
#

original_directory = os.getcwd()

# Determine number of runs.
# UGLY HACK
import deprecated_commands
output = deprecated_commands.getoutput('ls -1 %s | grep RUN' % project_directory)
elements = output.split()
nruns = len(elements)

for run in range(rank,nruns,size):
    run_directory = os.path.join(project_directory, 'RUN%d' % run)

    if not os.path.exists(run_directory):
        raise Exception("Expected RUN directory %s does not exist" % run_directory)

    # Determine number of clones.
    output = deprecated_commands.getoutput('ls -1 %s | grep state | grep xml' % run_directory)
    elements = output.split()
    nclones = len(elements)

    # Read System.
    system_xml_filename = os.path.join(run_directory, 'system.xml')
    system_xml = read_file(system_xml_filename)
    system = mm.XmlSerializer.deserialize(system_xml)

    for clone in range(nclones):
        clone_filename = os.path.join(run_directory, 'state%d.xml' % clone)
        if not os.path.exists(clone_filename):
            raise Exception("Expected CLONE file %s does not exist" % clone_filename)
    
        # Read Integrator.
        integrator_xml_filename = os.path.join(run_directory, 'integrator.xml')
        integrator_xml = read_file(integrator_xml_filename)
        integrator = mm.XmlSerializer.deserialize(integrator_xml)

        # Read state.
        state_xml = read_file(clone_filename)
        serialized_state = mm.XmlSerializer.deserialize(state_xml)
        
        # Create Context.
        context = mm.Context(system, integrator, platform)
        context.setPositions(serialized_state.getPositions())
        context.setVelocities(serialized_state.getVelocities())
        box_vectors = serialized_state.getPeriodicBoxVectors()
        context.setPeriodicBoxVectors(*box_vectors)

        # Check initial state.
        import numpy
        state = context.getState(getForces=True, getEnergy=True)
        #print "RUN%d CLONE%d: serialized State potential %s kinetic %s | computed potential from serialized positions %s" % (run, clone, str(serialized_state.getPotentialEnergy()), str(state.getKineticEnergy()), str(state.getPotentialEnergy()))

        # Check initial energy from serialized coordinates.
        if numpy.isnan(state.getPotentialEnergy() / units.kilocalories_per_mole):
            raise Exception("Initial energy is NaN after %d steps of dynamics starting from %s." % (nsteps, clone_filename))

        # TODO: Test forces.

        # Integrate.
        integrator.step(nsteps)

        # Check final state.
        import numpy
        state = context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True)
        #print "RUN%d CLONE%d: potential after dynamics is%s" % (run, clone, str(state.getPotentialEnergy()))
        if numpy.isnan(state.getPotentialEnergy() / units.kilocalories_per_mole):
            raise Exception("RUN%d CLONE%d: Final energy is NaN after %d steps of dynamics starting from %s." % (run, clone, nsteps, clone_filename))

        # Clean up.
        del context, state, integrator


