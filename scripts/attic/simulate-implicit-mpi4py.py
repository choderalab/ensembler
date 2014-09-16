# Subject models to implicit solvent simulation.
#
# John D. Chodera <choderaj@mskcc.org> - 17 Feb 2013
#
# PREREQUISITES
#
# * OpenMM
# http://simtk.org/home/openmm
#
# * mpi4py
# http://www.sagemath.org/doc/numerical_sage/mpi4py.html
#
# NOTE
# Must be launched with same number of CPU processes as GPUs, and processes must be assigned sequentially to each node (e.g. 0 1 2 to node 0, 3 4 5 to node 1, ...)

# PARAMETERS

import sys, traceback
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

process_only_these_templates = False
#process_only_these_templates = ['SRC_HUMAN_P12931_PK0_2H8H_A', 'HCK_HUMAN_P08631_PK0_1AD5_A', 'ABL1_HUMAN_P00519_PK0_2HZ0_A', 'ABL1_HUMAN_P00519_PK0_2HIW_A', 'CDK2_HUMAN_P24941_PK0_1KE7_A']

# OpenMM parameters

import simtk.openmm as mm
import simtk.unit as units

forcefields_to_use = ['amber99sbildn.xml', 'amber99_obc.xml'] # list of forcefields to use in parameterization

timestep = 2.0 * units.femtoseconds # timestep 
temperature = 300.0 * units.kelvin # simulation temperature 
collision_rate = 20.0 / units.picoseconds # Langevin collision rate
nsteps_per_iteration = 500 # number of timesteps per iteration
niterations = 100 # number of iterations
cutoff = None # nonbonded cutoff

minimization_tolerance = 10.0 * units.kilojoules_per_mole / units.nanometer
minimization_steps = 20

platform_name = 'CUDA'

kB = units.MOLAR_GAS_CONSTANT_R
kT = kB * temperature

pH = 8.0

verbose = False
write_trajectory = False

try:
        ngpus = int( sys.argv[ sys.argv.index('-gpupn') + 1 ] )
except ValueError:
        ngpus = 1 # number of GPUs per node

#
# GET ABSOLUTE PATHS
#

import os.path

# Input files.
targets_directory = os.path.abspath("targets") # target sequences for modeling
templates_directory = os.path.abspath("templates") # template structures for use in modeling
models_directory = os.path.abspath("models")

original_directory = os.getcwd()

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

import simtk.openmm.app as app
forcefield = app.ForceField(*forcefields_to_use)

#
# MULTIPROCESSING
#

def simulate_implicit(model_directory, variants=None, gpuid=0, rank=0):
    print "-------------------------------------------------------------------------"
    print "gpuid %d rank %d directory %s" % (gpuid, rank, model_directory)
    print "-------------------------------------------------------------------------"

    # Choose platform.
    platform = mm.Platform.getPlatformByName(platform_name)

    # Set GPU id.
    platform.setPropertyDefaultValue('CudaDeviceIndex', '%d' % gpuid)
    platform.setPropertyDefaultValue('OpenCLDeviceIndex', '%d' % gpuid)

    # Only simulate models that are unique following filtering by clustering.
    unique_by_clustering = os.path.exists(os.path.join(model_directory, 'unique_by_clustering'))
    if not unique_by_clustering: return

    os.chdir(model_directory)

    # Check to make sure the initial model file is present.
    model_filename = os.path.join(model_directory, 'model.pdb')
    if not os.path.exists(model_filename): return

    # Pass if this simulation has already been run.
    pdb_filename = os.path.join(model_directory, 'implicit-refined.pdb')
    if os.path.exists(pdb_filename): return

    if verbose: print "Reading model..."
    pdb = app.PDBFile(model_filename)

    # Add missing protons.
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield, pH=pH, variants=variants)
    topology = modeller.getTopology()
    positions = modeller.getPositions()

    if verbose: print "Constructing System object..."
    if cutoff is None:
        system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    else:
        system = forcefield.createSystem(topology, nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=cutoff, constraints=app.HBonds)
        
    if verbose: print "Creating Context..."
    integrator = mm.LangevinIntegrator(temperature, collision_rate, timestep)
    context = mm.Context(system, integrator, platform)
    context.setPositions(positions)

    if verbose: print "Minimizing structure..."
    mm.LocalEnergyMinimizer.minimize(context, minimization_tolerance, minimization_steps)

    if write_trajectory:
        # Open trajectory for writing.
        if verbose: print "Opening trajectory for writing..."
        trajectory_filename = os.path.join(model_directory, 'implicit-trajectory.pdb')
        trajectory_outfile = open(trajectory_filename, 'w')
        app.PDBFile.writeHeader(topology, file=trajectory_outfile)

    # Open energy trajectory for writing
    energy_filename = os.path.join(model_directory, 'implicit-energies.txt')
    energy_outfile = open(energy_filename, 'w')
    energy_outfile.write('# iteration | simulation time (ps) | potential_energy (kT) | kinetic_energy (kT) | ns per day\n')

    if verbose: print "Running dynamics..."
    import time
    initial_time = time.time()
    for iteration in range(niterations):
        # integrate dynamics
        integrator.step(nsteps_per_iteration)
        # get current state
        state = context.getState(getEnergy=True, getPositions=True)
        simulation_time = state.getTime()
        potential_energy = state.getPotentialEnergy()
        kinetic_energy = state.getKineticEnergy()
        final_time = time.time()
        elapsed_time = (final_time - initial_time) * units.seconds
        ns_per_day = (simulation_time / elapsed_time) / (units.nanoseconds / units.day)
        if verbose: print "  %8.1f ps : potential %8.3f kT | kinetic %8.3f kT | %.3f ns/day | %.3f s remain" % (simulation_time / units.picoseconds, potential_energy / kT, kinetic_energy / kT, ns_per_day, elapsed_time * (niterations-iteration-1) / (iteration+1) / units.seconds)

        # Check energies are still finite.
        import numpy
        if numpy.isnan(potential_energy/kT) or numpy.isnan(kinetic_energy/kT):
            raise Exception("Potential or kinetic energies are nan.")

        if write_trajectory:
            app.PDBFile.writeModel(topology, state.getPositions(), file=trajectory_outfile, modelIndex=iteration)
                
        # write data
        energy_outfile.write("  %8d %8.1f %8.3f %8.3f %.3f\n" % (iteration, simulation_time / units.picoseconds, potential_energy / kT, kinetic_energy / kT, ns_per_day))
        energy_outfile.flush()

    if write_trajectory:
        app.PDBFile.writeFooter(topology, file=trajectory_outfile)
        trajectory_outfile.close()

    energy_outfile.close()

    # Write final PDB file.
    pdb_outfile = open(pdb_filename, 'w')
    app.PDBFile.writeHeader(topology, file=pdb_outfile)
    app.PDBFile.writeFile(topology, state.getPositions(), file=pdb_outfile)
    app.PDBFile.writeFooter(topology, file=pdb_outfile)
    pdb_outfile.close()

    os.chdir(original_directory)    

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

# iterate over our subset
for target in targets:

    # Process only specified targets if directed.
    if process_only_these_targets and (target not in process_only_these_targets): continue

    target_directory = os.path.join(models_directory, target)
    if not os.path.exists(target_directory): continue

    #
    # DETERMINE PROTONATION STATE TO USE THROUGHOUT
    #
    
    # Determine highest-identity model.
    sequence_identities_filename = os.path.join(target_directory, 'sequence-identities.txt')
    infile = open(sequence_identities_filename, 'r')
    contents = infile.readline() # first line is highest sequence identity
    infile.close()
    [reference_template, reference_identity] = contents.split()
    if verbose: print "Using %s as highest identity model (%s%%)" % (reference_template, reference_identity)
    
    # Read PDB for reference model.
    reference_pdb_filename = os.path.join(target_directory, reference_template, 'model.pdb')
    reference_pdb = app.PDBFile(reference_pdb_filename)

    # Add missing protons.
    modeller = app.Modeller(reference_pdb.topology, reference_pdb.positions)
    variants = modeller.addHydrogens(forcefield, pH=pH)
    if verbose: 
        print "Reference variants extracted:"
        if variants != None:
            for (residue_index, residue) in enumerate(variants):
                if residue != None:
                    print "%8d %s" % (residue_index+1, residue)
            print ""
        else: print variants

    for template_index in range(rank, len(templates), size):
        template = templates[template_index]
        if process_only_these_templates and (template not in process_only_these_templates): continue

        #print "-------------------------------------------------------------------------"
        #print "Simulating %s => %s in implicit solvent for %.1f ps" % (target, template, niterations * nsteps_per_iteration * timestep / units.picoseconds)
        #print "-------------------------------------------------------------------------"
        
        model_directory = os.path.join(models_directory, target, template)
        if not os.path.exists(model_directory): continue

        try:
            simulate_implicit(model_directory, variants, gpuid, rank)
        except Exception as e:
            # Record rejected models.
            trbk = traceback.format_exc()
            reject_file_path = os.path.join(model_directory, 'implicit-rejected.txt')
            reject_file = open(reject_file_path, 'w')
            reject_file.write(trbk)
            reject_file.close()

