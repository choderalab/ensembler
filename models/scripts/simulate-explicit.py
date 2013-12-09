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
# Process only these targets, if specified.
# e.g. -targets '["SRC_HUMAN_PK0_P12931", "ABL1_HUMAN_PK0_P00519"]'
try:
    process_only_these_targets = sys.argv[ sys.argv.index('-targets') + 1 ]
else:
    process_only_these_targets = False

# OpenMM parameters

import simtk.openmm as mm
import simtk.unit as units
import simtk.openmm.app as app

forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml'] # list of forcefields to use in parameterization

timestep = 2.0 * units.femtoseconds # timestep 
temperature = 300.0 * units.kelvin # simulation temperature 
pressure = 1.0 * units.atmospheres # simulation pressure
collision_rate = 20.0 / units.picoseconds # Langevin collision rate
barostat_period = 50
nsteps_per_iteration = 500 # number of timesteps per iteration
niterations = 100 # number of iterations

nonbondedMethod = app.PME

minimization_tolerance = 10.0 * units.kilojoules_per_mole / units.nanometer
minimization_steps = 20

platform = mm.Platform.getPlatformByName('CUDA')

kB = units.MOLAR_GAS_CONSTANT_R
kT = kB * temperature

verbose = False
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

templates_index_filename = os.path.join(templates_directory, 'templates.txt')
infile = open(templates_index_filename, 'r')
templates = [ line.strip() for line in infile ]
infile.close()
print '%d template structures' % len(templates)

#
# LOAD FORCEFIELD
#

forcefield = app.ForceField(*forcefields_to_use)

#
# SOLVATION METHOD
#

def solvate_pdb(pdb, target_nwaters, model='tip3p'):
    """
    Solvate the contents of a PDB file, ensuring it has exactly 'target_nwaters' waters.

    ARGUMENTS
    
    pdb (simtk.openmm.app.PDBFile) - the PDB file to solvate
    nwaters (int) - number of waters to end up with

    OPTIONAL ARGUMENTS

    model (string) - solvent model to use (default: 'tip3p')

    RETURNS

    positions (list of list of simtk.unit.Quantity) - positions of particles
    topology (simtk.openmm.app.Topology) - topology object for solvated system

    ALGORITHM

    The system is initially solvated with a box of size 'boxsize_guess'.
    If the system has too few waters, the boxsize is scaled by boxsize_enlarge_factor
    Once a sufficient number of waters are present, the last few waters are deleted to ensure target_nwaters is achieved.

    TODO

    There is no error checking to be sure that waters are not initially present in the system or the initially-present molecules are not deleted.

    """

    natoms_per_solvent = 3

    # Count initial atoms.
    natoms_initial = len(pdb.positions)
    if verbose: print "System initially has %d atoms (0 waters)" % (natoms_initial)

    # Solvate with zero padding to determine min number of waters and minimal unit cell dimensions.
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield, model=model, padding=0.0*units.angstroms)    
    topology = modeller.getTopology()
    positions = modeller.getPositions()
    box_min = topology.getUnitCellDimensions()
    natoms_min = len(positions) # minimal number of atoms
    nwaters_min = (natoms_min - natoms_initial) / natoms_per_solvent # minimal number of waters
    volume_min = box_min[0] * box_min[1] * box_min[2]
    residues = [ r for r in topology.residues() ] # build a list of residues
    nresidues_min = len(residues) # number of residues
    print box_min
    if verbose: print "Minimally solvated system has %d atoms (%d waters)" % (natoms_min, nwaters_min)

    # If minimally solvated box has more than target number of waters, reject.
    if (nwaters_min > target_nwaters):
        raise Exception("Minimally solvated system has more waters than target number of waters.  This indicates protein structure too extended to be useful.  Rejecting.")

    # Increase the box size by 10% and resolvate.
    scale = 1.1
    modeller = app.Modeller(pdb.topology, pdb.positions)
    topology = modeller.getTopology()
    topology.setUnitCellDimensions(box_min * scale)
    modeller.addSolvent(forcefield, model=model)    
    positions = modeller.getPositions()
    box_enlarged = topology.getUnitCellDimensions()
    natoms_enlarged = len(positions) # minimal number of atoms
    nwaters_enlarged = (natoms_enlarged - natoms_initial) / natoms_per_solvent # minimal number of waters
    volume_enlarged = box_enlarged[0] * box_enlarged[1] * box_enlarged[2]
    density = (nwaters_enlarged - nwaters_min) / (volume_enlarged - volume_min)
    if verbose: print "Enlarged solvated system has %d atoms (%d waters) : density of %.3f waters / nm^3" % (natoms_enlarged, nwaters_enlarged, density / (1.0 / units.nanometer**3))
    
    # Aim for slightly waters more than target.
    over_target = False
    extra_nwaters = 100
    while not over_target:
        delta_volume = (target_nwaters + extra_nwaters - nwaters_min) / density 
        scale = ((volume_min + delta_volume) / volume_min)**(1.0/3.0)
        if verbose: print "Final target of %d waters, so attempting box size %s to achieve %d waters..." % (target_nwaters, str(box_min * scale), target_nwaters + extra_nwaters)
        delta_volume = (target_nwaters + extra_nwaters - nwaters_min) / density 
        modeller = app.Modeller(pdb.topology, pdb.positions)
        topology = modeller.getTopology()
        topology.setUnitCellDimensions(box_min * scale)
        modeller.addSolvent(forcefield, model=model)
        positions = modeller.getPositions()
        topology = modeller.getTopology()
        natoms = len(positions) # minimal number of atoms
        nwaters = (natoms - natoms_initial) / natoms_per_solvent # minimal number of waters
        if verbose: print "  actual %d waters" % nwaters
        if (nwaters > target_nwaters):
            over_target = True
        else:
            extra_nwaters += 100
        
    # Delete waters to achieve target.
    ndelete = nwaters - target_nwaters
    if (ndelete > 0):
        if verbose: print "Will delete %d waters..." % ndelete
        residues = [ r for r in topology.residues() ] # build a list of residues
        nresidues = len(residues)
        
        # Select a random subset to delete.
        import numpy.random
        indices = numpy.random.permutation(range(nresidues_min,nresidues))
        residues_to_delete = list()
        for index in indices[0:ndelete]:
            residues_to_delete.append(residues[index])

        modeller.delete(residues_to_delete)

        # Get topology and positions.
        topology = modeller.getTopology()
        positions = modeller.getPositions()

        # Count number of waters.
        natoms_final = len(positions)
        nwaters = (natoms_final - natoms_initial) / 3

    print "%d waters total." % nwaters

    if (nwaters != target_nwaters):
        raise Exception("Malfunction in solvate_pdb: nwaters = %d, target_nwaters = %d" % (nwaters, target_nwaters))

    return [positions, topology]

#
# SIMULATE MODELS
#

original_directory = os.getcwd()

for target in targets:
    
    # Process only specified targets if directed.
    if process_only_these_targets and (target not in process_only_these_targets): continue

    target_directory = os.path.join(models_directory, target)
    if not os.path.exists(target_directory): continue

    # Start a 'reject file'.
    #reject_filename = os.path.join(target_directory, 'reject-explicit.txt')
    #reject_file = open(reject_filename, 'w')
    
    # Determine number of waters to use.
    nwaters_filename = os.path.join(target_directory, 'nwaters-use.txt')
    infile = open(nwaters_filename, 'r')
    line = infile.readline()
    nwaters = int(line)
    infile.close()

    # Process all templates.
    for template in templates:

        print "-------------------------------------------------------------------------"
        print "Simulating %s => %s in explicit solvent for %.1f ps" % (target, template, niterations * nsteps_per_iteration * timestep / units.picoseconds)
        print "-------------------------------------------------------------------------"
        
        model_directory = os.path.join(models_directory, target, template)
        if not os.path.exists(model_directory): continue

        # Only simulate models that are unique following filtering by clustering.
        unique_by_clustering = os.path.exists(os.path.join(model_directory, 'unique_by_clustering'))
        if not unique_by_clustering: continue

        os.chdir(model_directory)

        model_filename = os.path.join(model_directory, 'implicit-refined.pdb')
        if not os.path.exists(model_filename): continue

        try:
            if verbose: print "Reading model..."
            pdb = app.PDBFile(model_filename)

            if verbose: print "Solvating model to achieve target of %d waters..." % nwaters
            [positions, topology] = solvate_pdb(pdb, nwaters)

            if verbose: print "Constructing System object..."
            system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, constraints=app.HBonds)
            if verbose: print "  system has %d atoms" % (system.getNumParticles())

            # Add barostat.
            if verbose: print "Adding barostat..."
            barostat = mm.MonteCarloBarostat(pressure, temperature, barostat_period)
            system.addForce(barostat)

            if verbose: print "Creating Context..."
            integrator = mm.LangevinIntegrator(temperature, collision_rate, timestep)
            context = mm.Context(system, integrator, platform)
            context.setPositions(positions)

            if verbose: print "Minimizing structure..."
            mm.LocalEnergyMinimizer.minimize(context, minimization_tolerance, minimization_steps)

            if write_trajectory:
                # Open trajectory for writing.
                if verbose: print "Opening trajectory for writing..."
                trajectory_filename = os.path.join(model_directory, 'explicit-trajectory.pdb')
                trajectory_outfile = open(trajectory_filename, 'w')
                app.PDBFile.writeHeader(pdb.topology, file=trajectory_outfile)

            # Open energy trajectory for writing
            energy_filename = os.path.join(model_directory, 'explicit-energies.txt')
            energy_outfile = open(energy_filename, 'w')
            energy_outfile.write('# iteration | simulation time (ps) | potential_energy (kT) | kinetic_energy (kT) | volume (nm^3) | ns per day\n')
        
            if verbose: print "Running dynamics..."
            context.setVelocitiesToTemperature(temperature)
            import time
            initial_time = time.time()
            for iteration in range(niterations):
                # integrate dynamics
                integrator.step(nsteps_per_iteration)
                # get current state
                state = context.getState(getEnergy=True)
                simulation_time = state.getTime()
                potential_energy = state.getPotentialEnergy()
                kinetic_energy = state.getKineticEnergy()
                final_time = time.time()
                elapsed_time = (final_time - initial_time) * units.seconds
                ns_per_day = (simulation_time / elapsed_time) / (units.nanoseconds / units.day)
                box_vectors = state.getPeriodicBoxVectors()
                volume_in_nm3 = (box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]) / (units.nanometers**3) # TODO: Use full determinant
                remaining_time = elapsed_time * (niterations-iteration-1) / (iteration+1)
                if verbose: print "  %8.1f ps : potential %8.3f kT | kinetic %8.3f kT | volume %.3f nm^3 | %.3f ns/day | %.3f s remain" % (simulation_time / units.picoseconds, potential_energy / kT, kinetic_energy / kT, volume_in_nm3, ns_per_day, remaining_time / units.seconds)
            
                if write_trajectory:
                    state = context.getState(getPositions=True)
                    app.PDBFile.writeModel(pdb.topology, state.getPositions(), file=trajectory_outfile, modelIndex=iteration)
            
                # write data
                energy_outfile.write("  %8d %8.1f %8.3f %8.3f %.3f %.3f\n" % (iteration, simulation_time / units.picoseconds, potential_energy / kT, kinetic_energy / kT, volume_in_nm3, ns_per_day))
                energy_outfile.flush()

            if write_trajectory:
                app.PDBFile.writeFooter(pdb.topology, file=trajectory_outfile)
                trajectory_outfile.close()

            energy_outfile.close()

            state = context.getState(getPositions=True, enforcePeriodicBox=True)            
            pdb_filename = os.path.join(model_directory, 'explicit-refined.pdb')
            pdb_outfile = open(pdb_filename, 'w')
            app.PDBFile.writeHeader(topology, file=pdb_outfile)
            app.PDBFile.writeFile(topology, state.getPositions(), file=pdb_outfile)
            app.PDBFile.writeFooter(topology, file=pdb_outfile)
            pdb_outfile.close()

            # Serialize system
            if verbose: print "Serializing system..."
            system_filename = os.path.join(model_directory, 'explicit-system.xml')
            system_file = open(system_filename, 'w')
            system_file.write(mm.XmlSerializer.serialize(system))
            system_file.close()

            # Serialize integrator.
            if verbose: print "Serializing integrator..."
            integrator_filename = os.path.join(model_directory, 'explicit-integrator.xml')
            integrator_file = open(integrator_filename, 'w')
            integrator_file.write(mm.XmlSerializer.serialize(integrator))
            integrator_file.close()            

            # Serialize state.
            if verbose: print "Serializing state..."
            state = context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
            state_filename = os.path.join(model_directory, 'explicit-state.xml')
            state_file = open(state_filename, 'w')
            state_file.write(mm.XmlSerializer.serialize(state))
            state_file.close()            

            os.chdir(original_directory)    

        except Exception as e:
            import traceback
            print "Exception occurred with template %s" % template
            print traceback.format_exc()
            print str(e)

            # Add to rejection file.
            #reject_file.write('%s : %s\n' % (template, str(e)))
            #reject_file.flush()
            

