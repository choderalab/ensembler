#this is a map-function

class MSMSeed:
    
    alignment =None 
    seq_id = None
    model_struct = None
    failure_message = None
    implicit_refined_struct = None
    implicit_energies = None
    nwaters = None
    explicit_refined_struct = None
    system = None
    integrator = None
    state = None
    refinedModel = None
    explicitEnergies = None
    nwaters_min = 0
    natoms_min = 0
    volume_min = 0
    box_min = None
    tmpname = None
    tgtname =None
    def __init__(self, tgtseq, tmplseq,tmplstruct):
        self.tgtseq = tgtseq 
        self.tmplseq = tmplseq
        #this is just a string, since modeller is going to write it to file first
        self.tmplstruct = tmplstruct





#this function is an internal function that retrieves a specific chain from the PDB
#it returns the pdb structure of the chain as a string and the sequence as a string
#the last parameter is (mostly) for NMR structures with multiple models
#this is an internal function that is used by the first stage, which builds the MSMSeed object
def _retrieve_chain(pdb_code, chain_code, model_id=0):
    import Bio.PDB as pdb
    import tempfile
    import os
    import StringIO
    import urllib2
    get_fasta_restful = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=fastachain&compression=NO&structureId="
    query_string = get_fasta_restful+pdb_code + '&chainId=' + chain_code
    response = urllib2.urlopen(query_string)
    fasta_result = response.read()
    os.chdir(tempfile.mkdtemp())
    pdb_fetcher = pdb.PDBList()
    pdb_filepath = pdb_fetcher.retrieve_pdb_file(pdb_code)
    parser = pdb.PDBParser()
    structure = parser.get_structure(pdb_code, pdb_filepath)
    chain_result = structure[model_id][chain_code]
    outval = StringIO.StringIO()
    io = pdb.PDBIO()
    io.set_structure(chain_result)
    io.save(outval)
    return (fasta_result, outval.getvalue())


#perform alignment with clustal
def align_with_clustal(MSMSeedInput):
    import os.path
    import os
    import subprocess
    import tempfile
    import modeller
    import modeller.automodel
    #this...for now
    home_path = os.environ.get('MSMSEED_HOME')
    if home_path == None:
        home_path = '/cbio/jclab/pgrinaway'
         #The args object should contain two elements- a sequences string in seg format and a structure
    sequences = MSMSeedInput.tgtseq + MSMSeedInput.tmplseq

    
    #using modeller here to help with writing alignments exactly as it likes
    modeller.log.none()
    env = modeller.environ()         
         
    #the exit codes are useful for filtering RDDs, to reduce memory usage and improve performance
    MSMSeedInput.exitcode = -1
         
         #look for environment variable for clustal binary, otherwise use SD cluster location
    clustal_binary = os.environ.get('CLUSTAL_BINARY_PATH')
    if clustal_binary == None:
        clustal_path = '/cbio/jclab/share/kinome/external-tools/clustal-omega'
        clustal_binary_name = 'clustalo-1.1.0-linux-64'
        clustal_binary = os.path.join(clustal_path, clustal_binary_name)
         
         #make a temporary directory and change to it
    tmpath = tempfile.mkdtemp()
    os.chdir(tmpath)
         
    #write the unaligned sequences to a file in our tempdir, read them w/ the alignment obj
    #and rewrite them in fASTA format
    unaln_file = tempfile.NamedTemporaryFile(dir=tmpath, delete=False)
    unaln_file.writelines(sequences)
    unaln_file.close()
    aln = modeller.alignment(env)
    aln.append(unaln_file.name)
    aln.write('unaligned.fasta', alignment_format='FASTA')
         
         
    #We'll call clustal with commands, input- unaligned.fasta, output via STDOUT
    cmd = [clustal_binary, '--infile=unaligned.fasta','--dealign','--outfmt=vienna','--force']
    alignment = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin = subprocess.PIPE)
    aln_out = alignment.communicate()[0]
         
         
    #get the names of everything and their sequences separated
    tgtname = aln_out.splitlines()[0].strip('>\n')
    target_sequence = aln_out.splitlines()[1].strip('\n')
    tmpname= aln_out.splitlines()[2].strip('>\n')
    template_sequence = aln_out.splitlines()[3].strip('\n')

    #build up alignment file content string, but don't save to file (can do that in build_model)
    aln_pir = "Target-template alignment by clustal omega\n"
    aln_pir += ">P1;%s\n" % tgtname
    aln_pir += "sequence:SRC:FIRST:@:LAST :@:::-1.00:-1.00\n"
    aln_pir += target_sequence + '*\n' 
    aln_pir += ">P1;%s\n" % tmpname
    aln_pir += "structureX:%s:FIRST:@:LAST : :undefined:undefined:-1.00:-1.00\n" % tmpname
    aln_pir += template_sequence + '*\n'
    MSMSeedInput.alignment = aln_pir
    MSMSeedInput.tmpname = tmpname
    MSMSeedInput.tgtname = tgtname
    return MSMSeedInput 

#Function to perform alignments and build initial model structure. Requires file writing as Modeller does
#not at present accept objects

def build_model(MSMSeedInput):
         import os.path
         import os
         import tempfile
         import modeller
         import modeller.automodel
         import shutil
         from simtk.openmm import app
         home_path = os.environ.get('MSMSEED_HOME')
         if home_path == None:
              home_path = '/cbio/jclab/pgrinaway'
        
        
         structure = MSMSeedInput.tmplstruct
         tmpname = MSMSeedInput.tmpname
         tgtname = MSMSeedInput.tgtname
         #set an exit code- right now, 0 for success, 1 for identified failure, -1 for un
         #identified failure
         MSMSeedInput.exitcode = -2
         
         #look for environment variable for clustal binary, otherwise use SD cluster location
         
         
         #make a temporary directory and change to it
         tmpath = tempfile.mkdtemp()
         os.chdir(tmpath)
         
         #make a directory for the structure file
         struct_path = tempfile.mkdtemp()
         
         
         
         #initialize modeller, get environ object
         modeller.log.none()
         env = modeller.environ()
         
         env.io.atom_files_directory = ['.',struct_path]
         
         #write the structure of the template to a file with the name of the template, and tell env where to find it
         struct_file = open(os.path.join(struct_path,MSMSeedInput.tmpname),'w')
         struct_file.writelines(structure)
         struct_file.close()
         
         #write the alignmed sequences to file, modeller likes this
         aligned_filename = 'aligned.pir'
  
         outfile = open('aligned.pir', 'w')
         outfile.write(MSMSeedInput.alignment)
         outfile.close()
         #set up the automodel class
    
         a = modeller.automodel.allhmodel(env,
                                         # file with template codes and target sequence
                                         alnfile  = aligned_filename,
                                         # PDB codes of the template
                                         knowns   = tmpname,
                                         # code of the target
                                         sequence = tgtname)
                                         
        #make the model
         a.make()
         #if everything worked, collect the data and return, exit code 0
         if a.outputs[0]['failure']==None:
             tgt_model = modeller.model(env, file = a.outputs[0]['name'])
             tgt_model.write(file =  tgtname+tmpname)
             #write this to a pdb object
             pdb = app.PDBFile(tgtname+tmpname)
             MSMSeedInput.model_struct = pdb
            # model = open(tgtname+tmpname,'r').readlines()
             MSMSeedInput.seq_id = tgt_model.seq_id
             MSMSeedInput.exitcode=0
             #go home and wrap it up
             os.chdir(home_path)
             shutil.rmtree(tmpath)
             #that's it
             return MSMSeedInput
         else:
             #go home and wrap it up
             os.chdir(home_path)
             shutil.rmtree(tmpath)
             # didn't work :(
             MSMSeedInput.seq_id = tgt_model.seq_id
             MSMSeedInput.failure_msg = a.outputs[0]['failure']
             return MSMSeedInput

 
def simulate_implicit(modelSeed):
    import simtk.openmm as mm
    import simtk.unit as units
    import simtk.openmm.app as app
    import StringIO
    from py4j import JavaGateway
    #get stuff out of args.. should be the same format as output of modeller step:
     #  (exitcode, model, aln_out, seqid) = args
    variants=None
    #Make sure the previous step didn't fail (nonzero exit code)
    if modelSeed.exitcode!=0:
        return (-1,"Previous Stage Failed")
    
    
    #set some constants
    forcefields_to_use = ['amber99sbildn.xml', 'amber99_obc.xml'] # list of forcefields to use in parameterization 
    timestep = 2.0 * units.femtoseconds # timestep 
    temperature = 300.0 * units.kelvin # simulation temperature 
    collision_rate = 20.0 / units.picoseconds # Langevin collision rate
    nsteps_per_iteration = 500 # number of timesteps per iteration
    niterations = 100 # number of iterations
    cutoff = None # nonbonded cutoff
    minimization_tolerance = 10.0 * units.kilojoules_per_mole / units.nanometers
    minimization_steps = 20
    platform_name = 'CUDA'
    kB = units.MOLAR_GAS_CONSTANT_R
    kT = kB * temperature
    pH = 8.0
    verbose = False
    gpuid = 0
    #write_trajectory is not currently in use, but in the future may be useful
    write_trajectory = False
    forcefield = app.ForceField(*forcefields_to_use)
    platform = mm.Platform.getPlatformByName(platform_name)
    platform.setPropertyDefaultValue('CudaDeviceIndex', '%d' % gpuid)
    platform.setPropertyDefaultValue('OpenCLDeviceIndex', '%d' % gpuid)

    
    #open the model file
    pdb = modelSeed.model_struct
    
    #get modeller object
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield, pH=pH, variants=variants)
    topology = modeller.getTopology()
    positions = modeller.getPositions()
    
    
    #create system
    if verbose: print "Constructing System object..."
    if cutoff is None:
        system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    else:
        system = forcefield.createSystem(topology, nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=cutoff, constraints=app.HBonds)
        
    #create context    
    if verbose: print "Creating Context..."
    integrator = mm.LangevinIntegrator(temperature, collision_rate, timestep)
    context = mm.Context(system, integrator, platform)
    context.setPositions(positions)
    
    #do energy minimization!
    if verbose: print "Minimizing structure..."
    mm.LocalEnergyMinimizer.minimize(context, minimization_tolerance, minimization_steps)
    
    
    #no option to write trajectory for now--- would be bad for network performance
    #potentially want to transmit energies back in real time for debugging--consider for future
    energy_outfile = StringIO.StringIO()
    energy_outfile.write('# iteration | simulation time (ps) | potential_energy (kT) | kinetic_energy (kT) | ns per day\n')
    if verbose: print "Running dynamics..."
    import time
    initial_time = time.time()
    
    
    
    #do simulation
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
                energy_outfile.write("  %8d %8.1f %8.3f %8.3f %.3f\n" % (iteration, simulation_time / units.picoseconds, potential_energy / kT, kinetic_energy / kT, ns_per_day))
                modelSeed.exitcode = -2
                modelSeed.failure_message = "Failed at implicit solvent stage due to nan"
                energy_outfile.flush()
                modelSeed.implicit_energies = energy_outfile.getvalue()
                return modelSeed
               
                


           
            #not currently an option, but maybe in future:
            #if write_trajectory:
             #   app.PDBFile.writeModel(topology, state.getPositions(), file=trajectory_outfile, modelIndex=iteration)
            
            # write data
            energy_outfile.write("  %8d %8.1f %8.3f %8.3f %.3f\n" % (iteration, simulation_time / units.picoseconds, potential_energy / kT, kinetic_energy / kT, ns_per_day))
            energy_outfile.flush()
            modelSeed.implicit_energies = energy_outfile.getvalue()
            modelSeed.implicit_refined_struct = (topology, state.getPositions())

    return modelSeed

def getwaters_model(implicit_refined_MSMseed):
    import simtk.openmm as mm
    import simtk.unit as units
    import simtk.openmm.app as app
    import StringIO
    verbose = False

    model = 'tip3p'
    (topology, positions) = implicit_refined_MSMseed.implicit_refined_struct
    natoms_per_solvent = 3
    forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml']
    forcefield = app.ForceField(*forcefields_to_use)
    # Count initial atoms.
    natoms_initial = len(positions)
    if verbose: print "System initially has %d atoms (0 waters)" % (natoms_initial)

    # Solvate with zero padding to determine min number of waters and minimal unit cell dimensions.
    modeller = app.Modeller(topology, positions)
    modeller.addSolvent(forcefield, model='tip3p', padding=0.0*units.angstroms)    
    topology = modeller.getTopology()
    positions = modeller.getPositions()
    box_min = topology.getUnitCellDimensions()
    natoms_min = len(positions) # minimal number of atoms
    nwaters_min = (natoms_min - natoms_initial) / natoms_per_solvent # minimal number of waters
    volume_min = box_min[0] * box_min[1] * box_min[2]
   # residues = [ r for r in topology.residues() ] # build a list of residues
    #nresidues_min = len(residues) # number of residues
    implicit_refined_MSMseed.nwaters_min = nwaters_min
    implicit_refined_MSMseed.natoms_min = natoms_min
    implicit_refined_MSMseed.volume_min = volume_min
    implicit_refined_MSMseed.box_min = box_min
    implicit_refined_MSMseed.natoms_initial = natoms_initial

    
    return implicit_refined_MSMseed

#this is a utility function for getting a list of nwaters
#can be done more idiomatically with a reduce() but this can be
#useful for debugging
def get_nwaters_list(nwaters_determined_MSMseed):
    return nwaters_determined_MSMseed.nwaters_min

#given a list of numbers, and a percentile, returns the number at the pctile% of the list
def get_pctile(nwaters_list, pctile):
    import numpy
    nwaters_array = numpy.array(nwaters_list)
    nwaters_array.sort()
    pctile_idx =  int((len(nwaters_array) - 1) * pctile)
    return nwaters_array[pctile_idx]


#this function solvates the model
def solvate_model(implicit_refined_MSMseed):
    import simtk.openmm as mm
    import simtk.unit as units
    import simtk.openmm.app as app
    model = 'tip3p'
    verbose = True
    natoms_per_solvent = 3
    forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml']
    forcefield = app.ForceField(*forcefields_to_use)
    #the value below is broadcasted from the Spark master
    target_nwaters = 6000
    if (implicit_refined_MSMseed.nwaters_min > target_nwaters):
        implicit_refined_MSMseed.errorcode = -3
        implicit_refined_MSMseed.failure_message = "Minimally solvated system has more than the target number of waters"
        return implicit_refined_MSMseed

    # Increase the box size by 10% and resolvate.
    scale = 1.1
    (topology, positions) = implicit_refined_MSMseed.implicit_refined_struct
    modeller = app.Modeller(topology,positions)
    topology = modeller.getTopology()
    topology.setUnitCellDimensions(implicit_refined_MSMseed.box_min * scale)
    modeller.addSolvent(forcefield, model=model)    
    positions = modeller.getPositions()
    box_enlarged = topology.getUnitCellDimensions()
    residues = [ r for r in topology.residues() ] 
    nresidues_min = len(residues)
    natoms_enlarged = len(positions) # minimal number of atoms
    nwaters_enlarged = (natoms_enlarged - implicit_refined_MSMseed.natoms_initial) / natoms_per_solvent # minimal number of waters
    volume_enlarged = box_enlarged[0] * box_enlarged[1] * box_enlarged[2]
    density = (nwaters_enlarged - implicit_refined_MSMseed.nwaters_min) / (volume_enlarged - implicit_refined_MSMseed.volume_min)
    if verbose: print "Enlarged solvated system has %d atoms (%d waters) : density of %.3f waters / nm^3" % (natoms_enlarged, nwaters_enlarged, density / (1.0 / units.nanometer**3))
    
    # Aim for slightly waters more than target.
    over_target = False
    extra_nwaters = 100
    while not over_target:
        delta_volume = (target_nwaters + extra_nwaters - implicit_refined_MSMseed.nwaters_min) / density 
        scale = ((implicit_refined_MSMseed.volume_min + delta_volume) / implicit_refined_MSMseed.volume_min)**(1.0/3.0)
        if verbose: print "Final target of %d waters, so attempting box size %s to achieve %d waters..." % (target_nwaters, str(implicit_refined_MSMseed.box_min * scale), target_nwaters + extra_nwaters)
        delta_volume = (target_nwaters + extra_nwaters - implicit_refined_MSMseed.nwaters_min) / density 
        modeller = app.Modeller(topology, positions)
        topology = modeller.getTopology()
        topology.setUnitCellDimensions(implicit_refined_MSMseed.box_min * scale)
        modeller.addSolvent(forcefield, model=model)
        positions = modeller.getPositions()
        topology = modeller.getTopology()
        natoms = len(positions) # minimal number of atoms
        nwaters = (natoms - implicit_refined_MSMseed.natoms_initial) / natoms_per_solvent # minimal number of waters
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
        nwaters = (natoms_final - implicit_refined_MSMseed.natoms_initial) / 3

    print "%d waters total." % nwaters

    if (nwaters != target_nwaters):
        raise Exception("Malfunction in solvate_pdb: nwaters = %d, target_nwaters = %d" % (nwaters, target_nwaters))
    implicit_refined_MSMseed.solvatedModel = (positions, topology)
    
    return implicit_refined_MSMseed

def simulate_explicit(solvatedMSMSeed):
    import simtk.openmm as mm
    import simtk.unit as units
    import simtk.openmm.app as app
    import StringIO
    timestep = 2.0 * units.femtoseconds # timestep 
    temperature = 300.0 * units.kelvin # simulation temperature 
    pressure = 1.0 * units.atmospheres # simulation pressure
    collision_rate = 20.0 / units.picoseconds # Langevin collision rate
    barostat_period = 50
    nsteps_per_iteration = 500 # number of timesteps per iteration
    niterations = 100 # number of iterations
    write_trajectory = False
    nonbondedMethod = app.PME

    minimization_tolerance = 10.0 * units.kilojoules_per_mole / units.nanometer
    minimization_steps = 20

    platform = mm.Platform.getPlatformByName('CUDA')

    kB = units.MOLAR_GAS_CONSTANT_R
    kT = kB * temperature
    verbose =True
    (positions, topology) = solvatedMSMSeed.solvatedModel
    nwaters = 2000
    forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml']
    forcefield = app.ForceField(*forcefields_to_use)

        
   
    try:
            if verbose: print "Reading model..."
            

            if verbose: print "Solvating model to achieve target of %d waters..." % nwaters
            

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
                if verbose: print "Not implemented"
                
               

            # Open energy trajectory for writing
           
            energy_outfile = StringIO.StringIO()
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
                    print "not implemented"
            
                # write data
                energy_outfile.write("  %8d %8.1f %8.3f %8.3f %.3f %.3f\n" % (iteration, simulation_time / units.picoseconds, potential_energy / kT, kinetic_energy / kT, volume_in_nm3, ns_per_day))
                energy_outfile.flush()

            if write_trajectory:
                print "not implemented"

            energies = energy_outfile.getValue()
            solvatedMSMSeed.explicitEnergies = energies

            state = context.getState(getPositions=True, enforcePeriodicBox=True)            
            
            solvatedMSMSeed.refinedModel=  (topology, state.getPositions())
           

            solvatedMSMSeed.system = system


            solvatedMSMSeed.integrator = integrator           

            state = context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
         
            solvatedMSMSeed.state = state

            return solvatedMSMSeed 

    except Exception as e:
            import traceback
            print "Exception occurred with template"
            print traceback.format_exc()
            print str(e)

                       
        
             
             
            
         
         