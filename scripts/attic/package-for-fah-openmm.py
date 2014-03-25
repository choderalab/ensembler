# Package OpenMM serialized systems in Folding@Home project directory structure.
#
# John D. Chodera <choderaj@mskcc.org> - 16 Mar 2013

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

# Verbose output
verbose = True

# Number of clones
nclones = 10

# If True, try to tgz the whole thing.
archive = False

#
# PACKAGE RUN
#

def generateRun(project_directory, source_directory, run, nclones, verbose=False):
    """
    Build Folding@Home RUN and CLONE subdirectories from (possibly compressed) OpenMM serialized XML files.

    ARGUMENTS

    project_directory (string) - base project directory to place RUN in
    source_directory (string) - source directory for OpenMM serialized XML files
    run (int) - run index
    nclones (int) - number of clones to generate

    """

    if verbose: print "Building RUN %d" % run
 
    try:
        import simtk.openmm
        import os, os.path, shutil
        
        # Determine directory and pathnames.
        rundir = os.path.join(project_directory, 'RUN%d' % run)
        template_filename = os.path.join(rundir, 'template.txt')
        seqid_filename = os.path.join(rundir, 'sequence-identity.txt')
        system_filename = os.path.join(rundir, 'system.xml')
        integrator_filename = os.path.join(rundir, 'integrator.xml')
        protein_structure_filename = os.path.join(rundir, 'protein.pdb')
        system_structure_filename = os.path.join(rundir, 'system.pdb')
        protein_structure_filename_source = os.path.join(source_directory, 'implicit-refined.pdb')
        system_structure_filename_source = os.path.join(source_directory, 'explicit-refined.pdb')

        # Return if this directory has already been set up.
        if os.path.exists(rundir): 
            if os.path.exists(template_filename) and os.path.exists(seqid_filename) and os.path.exists(system_filename) and os.path.exists(integrator_filename) and os.path.exists(protein_structure_filename) and os.path.exists(system_structure_filename): return
        else:
            # Construct run directory if it does not exist.
            os.makedirs(rundir)

        # Write template information.
        [filepath, template_name] = os.path.split(source_directory)
        outfile = open(template_filename, 'w')
        outfile.write(template_name + '\n')
        outfile.close()

        # Copy the protein and system structure pdbs
        shutil.copyfile(protein_structure_filename_source, protein_structure_filename)
        shutil.copyfile(system_structure_filename_source, system_structure_filename)

        # Read system, integrator, and state.
        def readFileContents(filename):
            import os.path
            fullpath = os.path.join(source_directory, filename)

            if os.path.exists(fullpath):
                infile = open(fullpath, 'r')
            elif os.path.exists(fullpath+'.gz'):
                import gzip
                infile = gzip.open(fullpath+'.gz', 'r')
            else:
                raise IOError('File %s not found' % filename)

            contents = infile.read()
            infile.close()
            return contents
        def writeFileContents(filename, contents):
            outfile = open(filename, 'w')
            outfile.write(contents)
            outfile.close()

        import simtk.openmm
        system     = simtk.openmm.XmlSerializer.deserialize(readFileContents('explicit-system.xml'))
        state      = simtk.openmm.XmlSerializer.deserialize(readFileContents('explicit-state.xml'))

        # Substitute default box vectors.
        box_vectors = state.getPeriodicBoxVectors()
        system.setDefaultPeriodicBoxVectors(*box_vectors)

        # Write sequence identity.
        contents = readFileContents(os.path.join(source_directory, 'sequence-identity.txt'))
        writeFileContents(seqid_filename, contents)

        # Integrator settings.
        import simtk.unit as units
        constraint_tolerance = 1.0e-5 
        timestep = 2.0 * units.femtoseconds
        collision_rate = 5.0 / units.picosecond
        temperature = 300.0 * units.kelvin

        # Create new integrator to use.
        integrator = simtk.openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        
        # TODO: Make sure MonteCarloBarostat temperature matches set temperature.

        # Serialize System.
        writeFileContents(system_filename, simtk.openmm.XmlSerializer.serialize(system))

        # Serialize Integrator
        writeFileContents(integrator_filename, simtk.openmm.XmlSerializer.serialize(integrator))

        # Create Context so we can randomize velocities.
        platform = simtk.openmm.Platform.getPlatformByName('Reference')
        context = simtk.openmm.Context(system, integrator, platform)
        context.setPositions(state.getPositions())
        context.setVelocities(state.getVelocities())
        box_vectors = state.getPeriodicBoxVectors()
        context.setPeriodicBoxVectors(*box_vectors)
            
        # Create clones with different random initial velocities.
        for clone_index in range(nclones):
            context.setVelocitiesToTemperature(temperature)
            state = context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)                                     
            state_filename = os.path.join(rundir, 'state%d.xml' % clone_index)
            writeFileContents(state_filename, simtk.openmm.XmlSerializer.serialize(state))

        # Clean up.
        del context, integrator, state, system

    except Exception as e:
        import traceback
        print traceback.format_exc()
        print str(e)    

    return


if __name__ == '__main__':

    #
    # GET ABSOLUTE PATHS
    #

    import os.path

    # Input files.
    targets_directory = os.path.abspath("targets") # target sequences for modeling
    templates_directory = os.path.abspath("templates") # template structures for use in modeling
    models_directory = os.path.abspath("models")

    # Output files.
    projects_directory = os.path.abspath("projects") # FAH projects directory

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
    # SET UP PROJECTS
    #

    import os
    if not os.path.exists(projects_directory):
        os.makedirs(projects_directory)

    for target in targets:
        
        # Process only specified targets if directed.
        if process_only_these_targets and (target not in process_only_these_targets): continue

        target_directory = os.path.join(models_directory, target)
        if not os.path.exists(target_directory): continue

        print "-------------------------------------------------------------------------"
        print "Building FAH OpenMM project for target %s" % (target)
        print "-------------------------------------------------------------------------"

        #
        # BUILD A LIST OF VALID TEMPLATES
        #

        # Process all templates.
        if verbose: print "Building list of valid templates..."
        import os.path
        valid_templates = list()
        for template in templates:
            # Check to make sure all files needed are present.
            is_valid = True
            filenames = ['explicit-system.xml', 'explicit-state.xml', 'explicit-integrator.xml']
            for filename in filenames:
                fullpath = os.path.join(target_directory, template, filename)
                if not (os.path.exists(fullpath) or os.path.exists(fullpath+'.gz')):
                    is_valid = False
            # Exclude those that are not unique by clustering.
            unique_by_clustering = os.path.exists(os.path.join(target_directory, template, 'unique_by_clustering'))
            if not unique_by_clustering:
                is_valid = False
            # TODO: Exclude if final potential energies from explicit solvent equilibration are too high.

            # Append if valid.
            if is_valid:
                valid_templates.append(template)

        nvalid = len(valid_templates)
        if verbose: print "%d valid unique initial starting conditions found" % nvalid

        #
        # SORT BY SEQUENCE IDENTITY
        #

        if verbose: print "Sorting templates in order of decreasing sequence identity..."
        import numpy
        sequence_identities = numpy.zeros([nvalid], numpy.float32)
        for (template_index, template) in enumerate(valid_templates):
            filename = os.path.join(target_directory, template, 'sequence-identity.txt')
            infile = open(filename, 'r')
            contents = infile.readline().strip()
            infile.close()
            sequence_identity = float(contents)
            sequence_identities[template_index] = sequence_identity
        sorted_indices = numpy.argsort(-sequence_identities)
        valid_templates = [ valid_templates[index] for index in sorted_indices ]
        if verbose: 
            print "Sorted"
            print sequence_identities[sorted_indices]

        #
        # CREATE PROJECT DIRECTORY
        # 

        project_directory = os.path.join(projects_directory, target)
        if not os.path.exists(projects_directory):
            os.makedirs(projects_directory)

        #
        # BUILD RUNS IN PARALLEL
        #

        if verbose: print "Building RUNs in parallel..."
        import multiprocessing
        pool = multiprocessing.Pool()    
        results = list()
        for (run_index, template) in enumerate(valid_templates):
            source_directory = os.path.join(target_directory, template)
            result = pool.apply_async(generateRun, args=(project_directory, source_directory, run_index, nclones, verbose))    
            results.append(result)
        pool.close()
        pool.join()

        # 
        # ARCHIVE
        #

        if archive:
            if verbose: print "Generating archive ..."
            import commands
            archive_filename = os.path.join(projects_directory, target + '.tgz')
            project_directory = os.path.join(projects_directory, target)
            commands.getoutput('tar zcf %s %s' % (archive_filename, project_directory))
    


