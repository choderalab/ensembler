import ensembler


def package_for_fah(process_only_these_targets=None, verbose=False, nclones=10, archive=True):
    '''Create the input files and directory structure necessary to start a Folding@Home project.

    MPI-enabled.

    Parameters
    ----------
    archive : Bool
        A .tgz compressed archive will be created for each individual RUN directory.
    '''
    import os
    import numpy
    import mpi4py.MPI
    import simtk.openmm as openmm
    import simtk.unit as unit
    comm = mpi4py.MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    models_dir = os.path.abspath('models')
    packaged_models_dir = os.path.abspath('packaged-models')
    projects_dir = os.path.join(packaged_models_dir, 'fah-projects')
    original_dir = os.getcwd()
    if rank == 0:
        if not os.path.exists(projects_dir):
            os.mkdir(projects_dir)
    comm.Barrier()

    targets, templates_resolved_seq, templates_full_seq = ensembler.core.get_targets_and_templates()
    templates = templates_resolved_seq

    def generateRun(project_dir, source_dir, run, nclones, verbose=False):
        """
        Build Folding@Home RUN and CLONE subdirectories from (possibly compressed) OpenMM serialized XML files.

        ARGUMENTS

        project_dir (string) - base project directory to place RUN in
        source_dir (string) - source directory for OpenMM serialized XML files
        run (int) - run index
        nclones (int) - number of clones to generate

        """

        if verbose: print "Building RUN %d" % run
     
        try:
            import os, shutil
            import gzip
            
            # Determine directory and pathnames.
            rundir = os.path.join(project_dir, 'RUN%d' % run)
            template_filename = os.path.join(rundir, 'template.txt')
            seqid_filename = os.path.join(rundir, 'sequence-identity.txt')
            system_filename = os.path.join(rundir, 'system.xml')
            integrator_filename = os.path.join(rundir, 'integrator.xml')
            protein_structure_filename = os.path.join(rundir, 'protein.pdb')
            system_structure_filename = os.path.join(rundir, 'system.pdb')
            protein_structure_filename_source = os.path.join(source_dir, 'implicit-refined.pdb')
            system_structure_filename_source = os.path.join(source_dir, 'explicit-refined.pdb')

            # Return if this directory has already been set up.
            if os.path.exists(rundir): 
                if os.path.exists(template_filename) and os.path.exists(seqid_filename) and os.path.exists(system_filename) and os.path.exists(integrator_filename) and os.path.exists(protein_structure_filename) and os.path.exists(system_structure_filename): return
            else:
                # Construct run directory if it does not exist.
                os.makedirs(rundir)

            # Write template information.
            [filepath, template_name] = os.path.split(source_dir)
            with open(template_filename, 'w') as outfile:
                outfile.write(template_name + '\n')

            # Copy the protein and system structure pdbs
            shutil.copyfile(protein_structure_filename_source, protein_structure_filename)
            shutil.copyfile(system_structure_filename_source, system_structure_filename)

            # Read system, integrator, and state.
            def readFileContents(filename):
                fullpath = os.path.join(source_dir, filename)

                if os.path.exists(fullpath):
                    infile = open(fullpath, 'r')
                elif os.path.exists(fullpath+'.gz'):
                    infile = gzip.open(fullpath+'.gz', 'r')
                else:
                    raise IOError('File %s not found' % filename)

                contents = infile.read()
                infile.close()
                return contents

            def writeFileContents(filename, contents):
                with open(filename, 'w') as outfile:
                    outfile.write(contents)

            system     = openmm.XmlSerializer.deserialize(readFileContents('explicit-system.xml'))
            state      = openmm.XmlSerializer.deserialize(readFileContents('explicit-state.xml'))

            # Substitute default box vectors.
            box_vectors = state.getPeriodicBoxVectors()
            system.setDefaultPeriodicBoxVectors(*box_vectors)

            # Write sequence identity.
            contents = readFileContents(os.path.join(source_dir, 'sequence-identity.txt'))
            writeFileContents(seqid_filename, contents)

            # Integrator settings.
            constraint_tolerance = 1.0e-5 
            timestep = 2.0 * unit.femtoseconds
            collision_rate = 1.0 / unit.picosecond
            temperature = 300.0 * unit.kelvin

            # Create new integrator to use.
            integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
            
            # TODO: Make sure MonteCarloBarostat temperature matches set temperature.

            # Serialize System.
            writeFileContents(system_filename, openmm.XmlSerializer.serialize(system))

            # Serialize Integrator
            writeFileContents(integrator_filename, openmm.XmlSerializer.serialize(integrator))

            # Create Context so we can randomize velocities.
            platform = openmm.Platform.getPlatformByName('Reference')
            context = openmm.Context(system, integrator, platform)
            context.setPositions(state.getPositions())
            context.setVelocities(state.getVelocities())
            box_vectors = state.getPeriodicBoxVectors()
            context.setPeriodicBoxVectors(*box_vectors)
                
            # Create clones with different random initial velocities.
            for clone_index in range(nclones):
                context.setVelocitiesToTemperature(temperature)
                state = context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)                                     
                state_filename = os.path.join(rundir, 'state%d.xml' % clone_index)
                writeFileContents(state_filename, openmm.XmlSerializer.serialize(state))

            # Clean up.
            del context, integrator, state, system

        except Exception as e:
            import traceback
            print traceback.format_exc()
            print str(e)    

        return


    def archiveRun(run_index, verbose=False):
        import subprocess
        os.chdir(project_dir)
        archive_filename = 'RUN%d.tgz' % run_index
        run_dir = 'RUN%d' % run_index
        subprocess.call(['tar', 'zcf', archive_filename, run_dir])
        os.chdir(original_dir)




    for target in targets:

        # Process only specified targets if directed.
        if process_only_these_targets and (target.id not in process_only_these_targets): continue

        models_target_dir = os.path.join(models_dir, target.id)
        if rank == 0:
            if not os.path.exists(models_target_dir): continue

        comm.Barrier()

        if rank == 0:
            print "-------------------------------------------------------------------------"
            print "Building FAH OpenMM project for target %s" % target.id
            print "-------------------------------------------------------------------------"

        # ========
        # Build a list of valid templates
        # ========

        # Process all templates.
        if verbose: print "Building list of valid templates..."
        valid_templates = list()
        for template in templates:
            # Check to make sure all files needed are present.
            is_valid = True
            filenames = ['explicit-system.xml', 'explicit-state.xml', 'explicit-integrator.xml']
            for filename in filenames:
                fullpath = os.path.join(models_target_dir, template.id, filename)
                if not (os.path.exists(fullpath) or os.path.exists(fullpath+'.gz')):
                    is_valid = False
            # Exclude those that are not unique by clustering.
            unique_by_clustering = os.path.exists(os.path.join(models_target_dir, template.id, 'unique_by_clustering'))
            if not unique_by_clustering:
                is_valid = False
            # TODO: Exclude if final potential energies from explicit solvent equilibration are too high.

            # Append if valid.
            if is_valid:
                valid_templates.append(template)

        nvalid = len(valid_templates)
        if verbose: print "%d valid unique initial starting conditions found" % nvalid

        # ========
        # Sort by sequence identity
        # ========

        if verbose: print "Sorting templates in order of decreasing sequence identity..."
        sequence_identities = numpy.zeros([nvalid], numpy.float32)
        for (template_index, template) in enumerate(valid_templates):
            filename = os.path.join(models_target_dir, template.id, 'sequence-identity.txt')
            with open(filename, 'r') as infile:
                contents = infile.readline().strip()
            sequence_identity = float(contents)
            sequence_identities[template_index] = sequence_identity
        sorted_indices = numpy.argsort(-sequence_identities)
        valid_templates = [ valid_templates[index] for index in sorted_indices ]
        if verbose: 
            print "Sorted"
            print sequence_identities[sorted_indices]

        # ========
        # Create project directory
        # ========

        project_dir = os.path.join(projects_dir, target.id)
        if rank == 0:
            if not os.path.exists(project_dir):
                os.makedirs(project_dir)

        comm.Barrier()

        # ========
        # Build runs in parallel
        # ========

        if verbose: print "Building RUNs in parallel..."
        for run_index in range(rank, len(valid_templates), size):
            print "-------------------------------------------------------------------------"
            print "Building RUN for template %s" % valid_templates[run_index].id
            print "-------------------------------------------------------------------------"

            source_dir = os.path.join(models_target_dir, valid_templates[run_index].id)
            generateRun(project_dir, source_dir, run_index, nclones, verbose)
            if archive:
                archiveRun(run_index, verbose)

        # TODO - get this working

        # if rank == 0:
        #
        #     # ========
        #     # Metadata
        #     # ========
        #
        #     import sys
        #     import yaml
        #     import ensembler.version
        #     import simtk.openmm.version
        #     datestamp = ensembler.core.get_utcnow_formatted()
        #
        #     meta_filepath = os.path.join(models_target_dir, 'meta.yaml')
        #     with open(meta_filepath) as meta_file:
        #         metadata = yaml.load(meta_file, Loader=ensembler.core.YamlLoader)
        #
        #     metadata['package_for_fah'] = {
        #         'target_id': target.id,
        #         'datestamp': datestamp,
        #         'python_version': sys.version.split('|')[0].strip(),
        #         'python_full_version': ensembler.core.literal_str(sys.version),
        #         'ensembler_version': ensembler.version.short_version,
        #         'ensembler_commit': ensembler.version.git_revision,
        #         'biopython_version': Bio.__version__,
        #         'openmm_version': simtk.openmm.version.short_version,
        #         'openmm_commit': simtk.openmm.version.git_revision
        #     }
        #
        #     meta_filepath = os.path.join(project_dir, 'meta.yaml')
        #     metadata = ensembler.core.ProjectMetadata(metadata)
        #     metadata.write(meta_filepath)

    comm.Barrier()
    if rank == 0:
        print 'Done.'


def package_for_transfer(process_only_these_targets=None):
    raise Exception, 'Not implemented yet.'
