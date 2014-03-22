def refine_implicitMD(openmm_platform='CUDA', gpupn=1, process_only_these_targets=None, process_only_these_templates=None, verbose=True, write_trajectory=False):
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
    targets = Bio.SeqIO.parse(targets_fasta_filename, 'fasta')
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
        print "-------------------------------------------------------------------------"
        print "gpuid %d rank %d directory %s" % (gpuid, rank, model_dir)
        print "-------------------------------------------------------------------------"

        # Choose platform.
        platform = openmm.Platform.getPlatformByName(openmm_platform)

        # Set GPU id.
        if openmm_platform == 'CUDA':
            platform.setPropertyDefaultValue('CudaDeviceIndex', '%d' % gpuid)
        elif openmm_platform == 'OpenCL':
            platform.setPropertyDefaultValue('OpenCLDeviceIndex', '%d' % gpuid)

        # Only simulate models that are unique following filtering by clustering.
        unique_by_clustering = os.path.exists(os.path.join(model_dir, 'unique_by_clustering'))
        if not unique_by_clustering: return

        os.chdir(model_dir)

        # Check to make sure the initial model file is present.
        model_filename = os.path.join(model_dir, 'model.pdb')
        if not os.path.exists(model_filename): return

        # Pass if this simulation has already been run.
        pdb_filename = os.path.join(model_dir, 'implicit-refined.pdb')
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

            try:
                simulate_implicitMD(model_dir, variants, gpuid, rank)
            except Exception:
                # Record rejected models.
                trbk = traceback.format_exc()
                reject_file_path = os.path.join(models_target_dir, 'implicit-rejected.txt')
                with open(reject_file_path, 'w') as reject_file:
                    reject_file.write(trbk)


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
    targets = Bio.SeqIO.parse(targets_fasta_filename, 'fasta')
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

            # Pass if this simulation has already been run.
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


def determine_nwaters(process_only_these_targets=None, process_only_these_templates=None, verbose=False):
    '''Determine distribution of nwaters and filter out models above the 68th percentile.
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
        targets = Bio.SeqIO.parse(targets_fasta_filename, 'fasta')
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

