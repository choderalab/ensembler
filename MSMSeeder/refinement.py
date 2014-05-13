def refine_implicitMD(openmm_platform='CUDA', gpupn=1, process_only_these_targets=None, process_only_these_templates=None, verbose=False, write_trajectory=False):
    '''Run MD refinement in implicit solvent.

    MPI-enabled.
    '''
    import os, traceback
    import Bio.SeqIO
    import simtk.openmm as openmm
    import simtk.unit as unit
    import simtk.openmm.app as app
    import mpi4py.MPI
    comm = mpi4py.MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    gpuid = (rank % gpupn)

    targets_dir = os.path.abspath("targets")
    templates_dir = os.path.abspath("templates")
    models_dir = os.path.abspath("models")
    original_dir = os.getcwd()

    targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
    targets = list( Bio.SeqIO.parse(targets_fasta_filename, 'fasta') )
    templates_fasta_filename = os.path.join(templates_dir, 'templates.fa')
    templates = list( Bio.SeqIO.parse(templates_fasta_filename, 'fasta') )

    # ========
    # Simulation parameters
    # ========

    forcefields_to_use = ['amber99sbildn.xml', 'amber99_obc.xml'] # list of forcefields to use in parameterization

    timestep = 2.0 * unit.femtoseconds # timestep 
    temperature = 300.0 * unit.kelvin # simulation temperature 
    collision_rate = 20.0 / unit.picoseconds # Langevin collision rate
    nsteps_per_iteration = 500 # number of timesteps per iteration
    niterations = 100 # number of iterations
    cutoff = None # nonbonded cutoff

    minimization_tolerance = 10.0 * unit.kilojoules_per_mole / unit.nanometer
    minimization_steps = 20

    kB = unit.MOLAR_GAS_CONSTANT_R
    kT = kB * temperature

    pH = 8.0

    forcefield = app.ForceField(*forcefields_to_use)

    def simulate_implicitMD(model_dir, variants=None, gpuid=0, rank=0, verbose=False):
        os.chdir(model_dir)

        # Choose platform.
        platform = openmm.Platform.getPlatformByName(openmm_platform)

        # Set GPU id.
        if openmm_platform == 'CUDA':
            platform.setPropertyDefaultValue('CudaDeviceIndex', '%d' % gpuid)
        elif openmm_platform == 'OpenCL':
            platform.setPropertyDefaultValue('OpenCLDeviceIndex', '%d' % gpuid)

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
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions)

        if verbose: print "Minimizing structure..."
        openmm.LocalEnergyMinimizer.minimize(context, minimization_tolerance, minimization_steps)

        if write_trajectory:
            # Open trajectory for writing.
            if verbose: print "Opening trajectory for writing..."
            trajectory_filename = os.path.join(model_dir, 'implicit-trajectory.pdb')
            trajectory_outfile = open(trajectory_filename, 'w')
            app.PDBFile.writeHeader(topology, file=trajectory_outfile)

        # Open energy trajectory for writing
        energy_filename = os.path.join(model_dir, 'implicit-energies.txt')
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
            elapsed_time = (final_time - initial_time) * unit.seconds
            ns_per_day = (simulation_time / elapsed_time) / (unit.nanoseconds / unit.day)
            if verbose: print "  %8.1f ps : potential %8.3f kT | kinetic %8.3f kT | %.3f ns/day | %.3f s remain" % (simulation_time / unit.picoseconds, potential_energy / kT, kinetic_energy / kT, ns_per_day, elapsed_time * (niterations-iteration-1) / (iteration+1) / unit.seconds)

            # Check energies are still finite.
            import numpy
            if numpy.isnan(potential_energy/kT) or numpy.isnan(kinetic_energy/kT):
                raise Exception("Potential or kinetic energies are nan.")

            if write_trajectory:
                app.PDBFile.writeModel(topology, state.getPositions(), file=trajectory_outfile, modelIndex=iteration)
                    
            # write data
            energy_outfile.write("  %8d %8.1f %8.3f %8.3f %.3f\n" % (iteration, simulation_time / unit.picoseconds, potential_energy / kT, kinetic_energy / kT, ns_per_day))
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

        os.chdir(original_dir)    



    for target in targets:
        if process_only_these_targets and (target.id not in process_only_these_targets): continue
        models_target_dir = os.path.join(models_dir, target.id)
        if not os.path.exists(models_target_dir): continue

        # ========
        # Determine protonation state to use throughout
        # ========
        
        # Determine highest-identity model.
        seqids_filename = os.path.join(models_target_dir, 'sequence-identities.txt')
        with open(seqids_filename, 'r') as seqids_file:
            contents = seqids_file.readline() # first line is highest sequence identity
        [reference_template, reference_identity] = contents.split()
        if verbose: print "Using %s as highest identity model (%s%%)" % (reference_template, reference_identity)
        
        # Read PDB for reference model.
        reference_pdb_filename = os.path.join(models_target_dir, reference_template, 'model.pdb')
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
            if process_only_these_templates and (template.id not in process_only_these_templates): continue

            print "-------------------------------------------------------------------------"
            print "Simulating %s => %s in implicit solvent for %.1f ps" % (target.id, template.id, niterations * nsteps_per_iteration * timestep / unit.picoseconds)
            print "-------------------------------------------------------------------------"
            
            model_dir = os.path.join(models_target_dir, template.id)
            if not os.path.exists(model_dir): continue

            # Only simulate models that are unique following filtering by clustering.
            unique_by_clustering = os.path.exists(os.path.join(model_dir, 'unique_by_clustering'))
            if not unique_by_clustering: continue

            # Check to make sure the initial model file is present.
            model_filename = os.path.join(model_dir, 'model.pdb')
            if not os.path.exists(model_filename): continue

            # Pass if this simulation has already been run.
            pdb_filename = os.path.join(model_dir, 'implicit-refined.pdb')
            if os.path.exists(pdb_filename): continue

            try:
                simulate_implicitMD(model_dir, variants, gpuid, rank)
            except Exception:
                # Record rejected models.
                trbk = traceback.format_exc()
                reject_file_path = os.path.join(models_target_dir, 'implicit-rejected.txt')
                with open(reject_file_path, 'w') as reject_file:
                    reject_file.write(trbk)

    comm.Barrier()
    if rank == 0:
        print 'Done.'


def solvate_models(process_only_these_targets=None, process_only_these_templates=None, verbose=False, write_trajectory=False):
    '''Solvate models which have been through MD refinement with implict solvent.

    MPI-enabled.
    '''
    import os
    import Bio.SeqIO
    import simtk.unit as unit
    import simtk.openmm.app as app
    import mpi4py.MPI
    comm = mpi4py.MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    targets_dir = os.path.abspath("targets")
    templates_dir = os.path.abspath("templates")
    models_dir = os.path.abspath("models")
    original_dir = os.getcwd()

    targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
    targets = list( Bio.SeqIO.parse(targets_fasta_filename, 'fasta') )
    templates_fasta_filename = os.path.join(templates_dir, 'templates.fa')
    templates = list( Bio.SeqIO.parse(templates_fasta_filename, 'fasta') )

    # OpenMM parameters

    forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml'] # list of forcefields to use in parameterization
    nparticles_per_water = 3 # number of particles per water molecule

    #box_width = 90.0 * unit.angstroms
    #boxsize = box_width * openmm.Vec3(1,1,1)
    padding = 10.0 * unit.angstroms

    forcefield = app.ForceField(*forcefields_to_use)

    for target in targets:
        
        # Process only specified targets if directed.
        if process_only_these_targets and (target.id not in process_only_these_targets): continue

        models_target_dir = os.path.join(models_dir, target.id)
        if not os.path.exists(models_target_dir): continue

        # Start a 'reject file'.
        reject_filename = os.path.join(models_target_dir, 'reject-solvation.txt')
        reject_file = open(reject_filename, 'w')

        # Process all templates.
        for template_index in range(rank, len(templates), size):
            template = templates[template_index]

            print "-------------------------------------------------------------------------"
            print "Solvating %s => %s in explicit solvent" % (target.id, template.id)
            print "-------------------------------------------------------------------------"
            
            model_dir = os.path.join(models_target_dir, template.id)
            if not os.path.exists(model_dir): continue

            model_filename = os.path.join(model_dir, 'implicit-refined.pdb')
            if not os.path.exists(model_filename): continue

            # Pass if solvation has already been run for this model.
            nwaters_filename = os.path.join(model_dir, 'nwaters.txt')
            if os.path.exists(nwaters_filename): continue

            os.chdir(model_dir)

            try:
                if verbose: print "Reading model..."
                pdb = app.PDBFile(model_filename)

                # Count initial atoms.
                natoms_initial = len(pdb.positions)

                if verbose: print "Solvating model..."
                modeller = app.Modeller(pdb.topology, pdb.positions)
                #modeller.addSolvent(forcefield, model='tip3p', boxSize=boxsize)
                modeller.addSolvent(forcefield, model='tip3p', padding=padding)
                #topology = modeller.getTopology()
                positions = modeller.getPositions()

                # Count final atoms.
                natoms_final = len(positions)
                nwaters = (natoms_final - natoms_initial) / nparticles_per_water
                if verbose: print "Solvated model contains %d waters" % nwaters

                # Record waters.
                with open(nwaters_filename, 'w') as nwaters_file:
                    nwaters_file.write('%d\n' % nwaters)

                os.chdir(original_dir)    

            except Exception as e:
                # Add to rejection file.
                reject_file.write('%s : %s\n' % (template.id, str(e)))
                reject_file.flush()

    comm.Barrier()
    if rank == 0:
        print 'Done.'


def determine_nwaters(process_only_these_targets=None, process_only_these_templates=None, verbose=False):
    '''Determine distribution of nwaters, and select the value at the 68th percentile.
    '''
    import os, numpy
    import Bio.SeqIO
    import mpi4py.MPI
    comm = mpi4py.MPI.COMM_WORLD
    rank = comm.rank

    # Run serially
    if rank == 0:
        targets_dir = os.path.abspath("targets")
        templates_dir = os.path.abspath("templates")
        models_dir = os.path.abspath("models")

        targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
        targets = list( Bio.SeqIO.parse(targets_fasta_filename, 'fasta') )
        templates_fasta_filename = os.path.join(templates_dir, 'templates.fa')
        templates = list( Bio.SeqIO.parse(templates_fasta_filename, 'fasta') )

        for target in targets:

            # Process only specified targets if directed.
            if process_only_these_targets and (target.id not in process_only_these_targets): continue

            models_target_dir = os.path.join(models_dir, target.id)
            if not os.path.exists(models_target_dir): continue

            if verbose: print "Determining number of waters in each system from target '%s'..." % target.id

            nwaters_list = []
            for template in templates:

                model_dir = os.path.join(models_target_dir, template.id)
                if not os.path.exists(model_dir): continue

                try:
                    nwaters_filename = os.path.join(model_dir, 'nwaters.txt')
                    with open(nwaters_filename, 'r') as nwaters_file:
                        firstline = nwaters_file.readline()
                    nwaters = int(firstline)
                    nwaters_list.append(nwaters)

                except Exception:
                    pass

            nwaters_array = numpy.array(nwaters_list)
            nwaters_array.sort()

            nwaters_list_filename = os.path.join(models_target_dir, 'nwaters-list.txt')
            with open(nwaters_list_filename, 'w') as nwaters_list_file:
                for nwaters in nwaters_array:
                    nwaters_list_file.write('%12d\n' % nwaters)

            # display statistics
            index68 = int((len(nwaters_array) - 1) * 0.68)
            index95 = int((len(nwaters_array) - 1) * 0.95)
            print "min = %d, max = %d, mean = %.1f, 68%% = %.0f, 95%% = %.0f\n" % (nwaters_array.min(), nwaters_array.max(), nwaters_array.mean(), nwaters_array[index68], nwaters_array[index95])

            filename = os.path.join(models_target_dir, 'nwaters-max.txt')
            with open(filename, 'w') as outfile:
                outfile.write('%d\n' % nwaters_array.max())

            # Use 68th percentile.
            filename = os.path.join(models_target_dir, 'nwaters-use.txt')
            with open(filename, 'w') as outfile:
                outfile.write('%d\n' % nwaters_array[index68])

    comm.Barrier()
    if rank == 0:
        print 'Done.'


# ========
# MD refinement with explicit solvent
# ========

def refine_explicitMD(openmm_platform='CUDA', gpupn=1, process_only_these_targets=None, process_only_these_templates=None, verbose=False, write_trajectory=False, careful_cleaning=True):
    '''Run MD refinement in explicit solvent.

    MPI-enabled.
    '''
    import os, traceback
    import Bio.SeqIO
    import simtk.openmm as openmm
    import simtk.unit as unit
    import simtk.openmm.app as app
    import mpi4py.MPI
    comm = mpi4py.MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    gpuid = (rank % gpupn)

    targets_dir = os.path.abspath("targets")
    templates_dir = os.path.abspath("templates")
    models_dir = os.path.abspath("models")
    original_dir = os.getcwd()

    targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
    targets = list( Bio.SeqIO.parse(targets_fasta_filename, 'fasta') )
    templates_fasta_filename = os.path.join(templates_dir, 'templates.fa')
    templates = list( Bio.SeqIO.parse(templates_fasta_filename, 'fasta') )

    # ========
    # Simulation parameters
    # ========

    forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml'] # list of forcefields to use in parameterization

    timestep = 2.0 * unit.femtoseconds # timestep 
    temperature = 300.0 * unit.kelvin # simulation temperature 
    pressure = 1.0 * unit.atmospheres # simulation pressure 
    collision_rate = 20.0 / unit.picoseconds # Langevin collision rate
    barostat_period = 50
    nsteps_per_iteration = 500 # number of timesteps per iteration
    niterations = 100 # number of iterations

    nonbondedMethod = app.PME

    minimization_tolerance = 10.0 * unit.kilojoules_per_mole / unit.nanometer
    minimization_steps = 20

    kB = unit.MOLAR_GAS_CONSTANT_R
    kT = kB * temperature

    forcefield = app.ForceField(*forcefields_to_use)


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
        modeller.addSolvent(forcefield, model=model, padding=0.0*unit.angstroms)    
        topology = modeller.getTopology()
        positions = modeller.getPositions()
        box_min = topology.getUnitCellDimensions()
        natoms_min = len(positions) # minimal number of atoms
        nwaters_min = (natoms_min - natoms_initial) / natoms_per_solvent # minimal number of waters
        volume_min = box_min[0] * box_min[1] * box_min[2]
        residues = [ r for r in topology.residues() ] # build a list of residues
        nresidues_min = len(residues) # number of residues
        if verbose: print "Minimally solvated system has %d atoms (%d waters)" % (natoms_min, nwaters_min)

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
        if verbose: print "Enlarged solvated system has %d atoms (%d waters) : density of %.3f waters / nm^3" % (natoms_enlarged, nwaters_enlarged, density / (1.0 / unit.nanometer**3))

        # Aim for slightly more waters than target.
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

        if (nwaters != target_nwaters):
            raise Exception("Malfunction in solvate_pdb: nwaters = %d, target_nwaters = %d" % (nwaters, target_nwaters))

        return [positions, topology]


    def simulate_explicitMD(model_dir, gpuid=0, rank=0, verbose=False):
        import gzip
        os.chdir(model_dir)

        # Choose platform.
        platform = openmm.Platform.getPlatformByName(openmm_platform)

        # Set GPU id.
        if openmm_platform == 'CUDA':
            platform.setPropertyDefaultValue('CudaDeviceIndex', '%d' % gpuid)
        elif openmm_platform == 'OpenCL':
            platform.setPropertyDefaultValue('OpenCLDeviceIndex', '%d' % gpuid)

        if verbose: print "Constructing System object..."
        system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, constraints=app.HBonds)
        if verbose: print "  system has %d atoms" % (system.getNumParticles())

        # Add barostat.
        if verbose: print "Adding barostat..."
        barostat = openmm.MonteCarloBarostat(pressure, temperature, barostat_period)
        system.addForce(barostat)

        if verbose: print "Creating Context..."
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions)

        if verbose: print "Minimizing structure..."
        openmm.LocalEnergyMinimizer.minimize(context, minimization_tolerance, minimization_steps)

        if write_trajectory:
            # Open trajectory for writing.
            if verbose: print "Opening trajectory for writing..."
            trajectory_filename = os.path.join(model_dir, 'explicit-trajectory.pdb')
            trajectory_outfile = open(trajectory_filename, 'w')
            app.PDBFile.writeHeader(pdb.topology, file=trajectory_outfile)

        # Open energy trajectory for writing
        energy_filename = os.path.join(model_dir, 'explicit-energies.txt')
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
            elapsed_time = (final_time - initial_time) * unit.seconds
            ns_per_day = (simulation_time / elapsed_time) / (unit.nanoseconds / unit.day)
            box_vectors = state.getPeriodicBoxVectors()
            volume_in_nm3 = (box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]) / (unit.nanometers**3) # TODO: Use full determinant
            remaining_time = elapsed_time * (niterations-iteration-1) / (iteration+1)
            if verbose: print "  %8.1f ps : potential %8.3f kT | kinetic %8.3f kT | volume %.3f nm^3 | %.3f ns/day | %.3f s remain" % (simulation_time / unit.picoseconds, potential_energy / kT, kinetic_energy / kT, volume_in_nm3, ns_per_day, remaining_time / unit.seconds)

            if write_trajectory:
                state = context.getState(getPositions=True)
                app.PDBFile.writeModel(pdb.topology, state.getPositions(), file=trajectory_outfile, modelIndex=iteration)

            # write data
            energy_outfile.write("  %8d %8.1f %8.3f %8.3f %.3f %.3f\n" % (iteration, simulation_time / unit.picoseconds, potential_energy / kT, kinetic_energy / kT, volume_in_nm3, ns_per_day))
            energy_outfile.flush()

        if write_trajectory:
            app.PDBFile.writeFooter(pdb.topology, file=trajectory_outfile)
            trajectory_outfile.close()

        energy_outfile.close()

        state = context.getState(getPositions=True, enforcePeriodicBox=True)            
        with open(pdb_filename, 'w') as pdb_outfile:
            app.PDBFile.writeHeader(topology, file=pdb_outfile)
            app.PDBFile.writeFile(topology, state.getPositions(), file=pdb_outfile)
            app.PDBFile.writeFooter(topology, file=pdb_outfile)

        # Serialize system
        if verbose: print "Serializing system..."
        with gzip.open(system_filename+'.gz', 'w') as system_file:
            system_file.write(openmm.XmlSerializer.serialize(system))

        # Serialize integrator.
        if verbose: print "Serializing integrator..."
        with gzip.open(integrator_filename+'.gz', 'w') as integrator_file:
            integrator_file.write(openmm.XmlSerializer.serialize(integrator))

        # Serialize state.
        if verbose: print "Serializing state..."
        state = context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
        with gzip.open(state_filename+'.gz', 'w') as state_file:
            state_file.write(openmm.XmlSerializer.serialize(state))

        os.chdir(original_dir)    


    for target in targets:
        if process_only_these_targets and (target.id not in process_only_these_targets): continue
        models_target_dir = os.path.join(models_dir, target.id)
        if not os.path.exists(models_target_dir): continue

        # Determine number of waters to use.
        nwaters_filename = os.path.join(models_target_dir, 'nwaters-use.txt')
        with open(nwaters_filename, 'r') as infile:
            line = infile.readline()
        nwaters = int(line)

        for template_index in range(rank, len(templates), size):
            template = templates[template_index]

            model_dir = os.path.join(models_target_dir, template.id)
            if not os.path.exists(model_dir): continue

            # Only simulate models that are unique following filtering by clustering.
            unique_by_clustering = os.path.exists(os.path.join(model_dir, 'unique_by_clustering'))
            if not unique_by_clustering: continue

            # Check to make sure the initial model file is present.
            model_filename = os.path.join(model_dir, 'implicit-refined.pdb')
            if not os.path.exists(model_filename): continue

            # Check if explicit solvent results are already available and usable.
            pdb_filename = os.path.join(model_dir, 'explicit-refined.pdb')
            system_filename = os.path.join(model_dir, 'explicit-system.xml')
            integrator_filename = os.path.join(model_dir, 'explicit-integrator.xml')
            state_filename = os.path.join(model_dir, 'explicit-state.xml')
            if os.path.exists(pdb_filename) and (os.path.exists(system_filename) or os.path.exists(system_filename+'.gz')) and (os.path.exists(integrator_filename) or os.path.exists(integrator_filename+'.gz')) and (os.path.exists(state_filename) or os.path.exists(state_filename+'.gz')): 
                # If not using 'careful mode', just continue.
                if not careful_cleaning:
                    continue
                    
                # Check if we can deserialize explicit solvent files.
                try:
                    system     = openmm.XmlSerializer.deserialize(readFileContents(system_filename))
                    state      = openmm.XmlSerializer.deserialize(readFileContents(state_filename))
                    # Serialized objects are OK---skip.
                    continue
                except Exception as e:
                    print e
                    # Attempt to delete the problematic files.
                    for filename in [system_filename, state_filename]:
                        try:
                            os.remove(filename)
                        except:
                            print e
                    # Now try to solvate and simulate model.
                    pass

            print "-------------------------------------------------------------------------"
            print "Simulating %s => %s in explicit solvent for %.1f ps" % (target.id, template.id, niterations * nsteps_per_iteration * timestep / unit.picoseconds)
            print "-------------------------------------------------------------------------"

            try:
                if verbose: print "Reading model..."
                pdb = app.PDBFile(model_filename)

                if verbose: print "Solvating model to achieve target of %d waters..." % nwaters
                [positions, topology] = solvate_pdb(pdb, nwaters)

                simulate_explicitMD(model_dir, gpuid, rank)

            except Exception:
                # Record rejected models.
                trbk = traceback.format_exc()
                reject_file_path = os.path.join(models_target_dir, 'explicit-rejected.txt')
                with open(reject_file_path, 'w') as reject_file:
                    reject_file.write(trbk)

    comm.Barrier()
    if rank == 0:
        print 'Done.'


def readFileContents(filename):
    import os.path

    if os.path.exists(filename):
        infile = open(filename, 'r')
    elif os.path.exists(filename+'.gz'):
        import gzip
        infile = gzip.open(filename+'.gz', 'r')
    else:
        raise IOError('File %s not found' % filename)

    contents = infile.read()
    infile.close()
    return contents

