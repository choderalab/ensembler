__author__ = 'Patrick B. Grinaway'



class MSMSeed(object):
    """
    Container class for information used in creating starting structures with MSMseeder

    Parameters
    ----------

    target_sequence : String
        FASTA format string with the target amino acid sequence and name of target after >
    template_sequence : String
        FASTA format string with the template amino acid sequence and name of template after >
    template_structure : String
        contains the string representation of the template pdb structure
    blast_eval : float, optional
        contains the e-value obtained from blast when aligning the template and target. Default 0


    Attributes
    ----------
    blast_eval : float
        contains the e-value obtained via blast search.
    rmsd_to_reference : float
        contains the RMSD to a reference structure
    target_sequence : String
        containing the sequence of this model's target in fasta
    target_ID : String
        containing the name of this model's target
    template_ID : String
        containing the name of this model's template
    template_sequence : String
        containing the sequence of this model's template in fasta
    template_structure : simtk.openmm.app.PDBFile
        object containing the structure of this model's template
    alignment : String
        containing the PIR format (MODELLER style) alignment of the target and template
    target_model : simtk.openmm.app.PDBFile
        object containing the output of MODELLER
    target_restraints: String
        containing .rsr file output of MODELLER
    implicit_refined_model : MDSys
        object containing an openmm topology and positions of an implicitly-refined model
    error_message : String
        containing a brief explanation of the error condition on this model, if any
    error_state : int
        the stage at which the error occurred
    nwaters : int
        the number of waters used to solvate the model
    solvated_model : MDSys
        object containing an openmm topology and positions of the solvated model
    explicit_refined_pdb : String
        gzipped PDB file contents of explicitly-refined model
    explicit_refined_state : String
        gzipped xml serialized form of explicit-refined model state
    explicit_refined_system : String
        gzipped xml serialized form of explicit-refined model system
    explicit_refined_integrator : String
        gzipped xml serialized form of explicit-refined model integrator





    """

    def __init__(self, target_sequence, template_sequence, template_structure, blast_eval=0):
        import StringIO
        import Bio.SeqIO

        target_seq_record = list()
        template_seq_record = list()

        target_sequence_stringio = StringIO.StringIO(target_sequence)
        template_sequence_stringio = StringIO.StringIO(template_sequence)

        target_sequence_stringio.seek(0)
        target_sequence_stringio.seek(0)

        target_sequence_parser = Bio.SeqIO.parse(target_sequence_stringio,'fasta')
        template_sequence_parser = Bio.SeqIO.parse(template_sequence_stringio,'fasta')

        for record in target_sequence_parser:
            self._target_sequence= record.seq
            self._target_id = record.id
        for record in template_sequence_parser:
            self._template_sequence = record.seq
            self._template_id = record.id

        self._template_structure = template_structure
        self._error_state = 0
        self._nwaters = 0
        self._blast_eval = blast_eval


    @property
    def blast_eval(self):
        return self._blast_eval

    @property
    def rmsd_to_reference(self):
        return self._rmsd_to_reference
    @rmsd_to_reference.setter
    def rmsd_to_reference(self, rmsd):
        self._rmsd_to_reference = rmsd

    @property
    def target_sequence(self):
        return self._target_sequence
    @property
    def template_sequence(self):
        return self._template_sequence
    @property
    def template_structure(self):
        return self._template_structure

    @property
    def alignment(self):
        return self._alignment
    @alignment.setter
    def alignment(self, new_alignment):
        self._alignment = new_alignment

    @property
    def target_model(self):
        return self._target_model
    @target_model.setter
    def target_model(self, new_target_model):
        self._target_model = new_target_model

    @property
    def target_restraints(self):
        return self._target_restraints
    @target_restraints.setter
    def target_restraints(self,new_target_restraints):
        self._target_restraints = new_target_restraints

    @property
    def target_id(self):
        return self._target_id
    @property
    def template_id(self):
        return self._template_id

    @property
    def implicit_refined_model(self):
        return self._implicit_refined_model
    @implicit_refined_model.setter
    def implicit_refined_model(self, new_model):
        self._implicit_refined_model = new_model

    @property
    def error_state(self):
        return self._error_state
    @error_state.setter
    def error_state(self, error):
        self._error_state = error

    @property
    def error_message(self):
        return self._error_message
    @error_message.setter
    def error_message(self, message):
        self._error_message = message

    @property
    def nwaters(self):
        return self._nwaters
    @nwaters.setter
    def nwaters(self,num_waters):
        self._nwaters = num_waters


    @property
    def solvated_model(self):
        return self._solvated_model
    @solvated_model.setter
    def solvated_model(self, new_solvated_model):
        self._solvated_model = new_solvated_model

    @property
    def explicit_refined_pdb(self):
        return self._explicit_refined_pdb
    @explicit_refined_pdb.setter
    def explicit_refined_pdb(self, new_pdb):
        self._explicit_refined_pdb = new_pdb

    @property
    def explicit_refined_state(self):
        return self._explicit_refined_state
    @explicit_refined_state.setter
    def explicit_refined_state(self,new_state):
        self._explicit_refined_state = new_state

    @property
    def explicit_refined_integrator(self):
        return self._explicit_refined_integrator
    @explicit_refined_integrator.setter
    def explicit_refined_integrator(self, new_integrator):
        self._explicit_refined_integrator = new_integrator

    @property
    def explicit_refined_system(self):
        return self._explicit_refined_system
    @explicit_refined_system.setter
    def explicit_refined_system(self, new_system):
        self._explicit_refined_system = new_system



class MDSys(object):
    def __init__(self, topology, positions):
        self._topology = topology
        self._positions = positions

    @property
    def topology(self):
        return self._topology
    @property
    def positions(self):
        return self._positions


def blast_pdb_local(fasta_string, num_hits=1000):
    import subprocess
    import os
    import shlex
    blast_data = os.getenv("DATA_HOME")
    blast_query = 'blastp -db %s/pdbaa -max_target_seqs %d -outfmt' % (blast_data, num_hits)
    out_fmt = '7 qseqid sseqid evalue bitscore'
    blast_cmd = shlex.split(blast_query)
    blast_cmd.append(out_fmt)
    p = subprocess.Popen(blast_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    blast_aln, error = p.communicate(input=fasta_string)
    msmseeds = []
    for result in blast_aln:
        res_data = result.split("\t")
        e_value = float(res_data[2])
        template_chain_code =  "_".join(res_data[1].split("|")[3:])
        template_fasta, template_pdb = _retrieve_chain(template_chain_code)
        msmseeds.append(MSMSeed(fasta_string, template_fasta, template_pdb, e_value))
    return msmseeds







def _retrieve_chain(pdb_code_input, model_id=0):
    import Bio.PDB as pdb
    import Bio.Seq
    import tempfile
    import os
    import StringIO
    import urllib2
    import shutil
    import simtk.openmm.app as app
    pdb_code, chain_code = pdb_code_input.split("_")
    temp_dir = tempfile.mkdtemp()
    os.chdir(temp_dir)
    pdb_fetcher = pdb.PDBList()
    pdb_filepath = pdb_fetcher.retrieve_pdb_file(pdb_code)
    parser = pdb.PDBParser()
    structure = parser.get_structure(pdb_code, pdb_filepath)
    chain_result = structure[model_id][chain_code]
    pp = pdb.PPBuilder()
    seq = Bio.Seq.Seq('')
    for peptide in pp.build_peptides(chain_result):
        seq+=peptide.get_sequence()
    fasta_result = ">%s_%s\n" % (pdb_code, chain_code) + seq
    outval = StringIO.StringIO()
    io = pdb.PDBIO()
    io.set_structure(chain_result)
    io.save(outval)
    outval.seek(0)
    shutil.rmtree(temp_dir)
    return fasta_result, app.PDBFile(outval)

def align_template_to_reference(msmseed, ref_msmseed):
    import modeller
    import tempfile
    import shutil
    temp_dir = tempfile.mkdtemp()
    try:
        os.chdir(temp_dir)
        alignment_file = open('aln_tmp.pir','w')
        aln = _PIR_alignment(ref_msmseed.template_sequence, ref_msmseed.template_id, msmseed.template_sequence, msmseed.template_id)
        alignment_file.writelines(aln)
        alignment_file.close()
        template_file = open(msmseed.template_id + '.pdb','w')
        template_pdb = msmseed.template_structure
        template_pdb.writeFile(template_pdb.topology, template_pdb.positions, template_file)
        template_file.close()
        ref_pdb = ref_msmseed.template_structure
        ref_file = open(ref_msmseed.template_id + '.pdb', 'w')
        ref_pdb.writeFile(ref_pdb.topology, ref_pdb.positions, ref_file)
        ref_file.close()
        modeller.log.none()
        env = modeller.environ()
        env.io.atom_files_directory = temp_dir
        aln = modeller.alignment(env, file='aln_tmp.pir', align_codes=(ref_msmseed.template_id, msmseed.template_id))
        mdl  = modeller.model(env, file=ref_msmseed.template_id + '.pdb')
        mdl2 = modeller.model(env, file=msmseed.template_id+'.pdb')
        mdl.pick_atoms(aln, pick_atoms_set=1, atom_types='CA')
        x = mdl.superpose(mdl2, aln)
        return x, mdl
    finally:
        shutil.rmtree(temp_dir)
    return "=("



def blast_pdb(target_sequence, num_hits=1000):
    """
    Query the PDB using NCBI blast and return MSMSeeds initialized with the results

    Parameters
    ----------
    target_sequence : String
        The sequence of the target to use to query blast
    num_hits : int, optional
        The maximum number of hits returned by BLAST. Default: 1000

    Returns
    -------
    msmseeds : list of MSMSeed objects
        A list of MSMSeed objects initialized with a target sequence, template sequence, template structure,
        and BLAST e-value. Can be readily parallelized in Spark.
    """
    from Bio.Blast import NCBIWWW, NCBIXML
    result_handle = NCBIWWW.qblast("blastp", "pdb", target_sequence, hitlist_size=num_hits)
    blast_record = NCBIXML.read(result_handle)
    alignments = blast_record.alignments
    msmseeds = []
    for alignment in alignments:
        e_val = alignment.hsps[0].expect
        template_fasta, template_structure = _retrieve_chain(alignment.accession)
        msmseeds.append(MSMSeed(target_sequence,template_fasta, template_structure, e_val))
    return msmseeds


def _extract_seq(pdb):
    import mdtraj as md
    mdtop = md.Topology.from_openmm(pdb.topology)
    resilist = []
    for residue in mdtop.residues:
        resilist.append(residue)
    return resilist

#probably no longer needed
def _correct_template_fasta(template_fasta):
    #sorry
    return "".join(template_fasta.splitlines()[0].split("|")[0].replace(":","_")+"\n"+"".join(template_fasta.splitlines()[1:]))





def _PIR_alignment(target_sequence, target_id, template_sequence, template_id):
    import Bio.SubsMat
    import Bio.pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    import Bio.SeqIO

    matrix = matlist.gonnet
    gap_open = -10
    gap_extend = -0.5
    aln = Bio.pairwise2.align.globalds(target_sequence, template_sequence, matrix, gap_open, gap_extend)

    #put together PIR file
    contents = "Target-template alignment by clustal omega\n"
    contents += ">P1;%s\n" % target_id
    contents += "sequence:%s:FIRST:@:LAST :@:::-1.00:-1.00\n" % target_id
    contents += aln[0][0] + '*\n'
    contents += ">P1;%s\n" % template_id
    contents += "structureX:%s:FIRST:@:LAST : :undefined:undefined:-1.00:-1.00\n" % template_id
    contents += aln[0][1] + '*\n'
    return contents



def target_template_alignment(msmseed):
    """
    Use Biopython pairwise2 to generate a sequence alignment in PIR format so that MODELLER can be used.
    Puts alignment in MSMSeed. The functionality for this was separated into another internal function,
    _PIR_alignment().

    Parameters
    ----------

    msmseed : MSMSeed
        object containing the sequences to be aligned

    Returns
    -------
        msmseed : MSMSeed
            Object containing the alignment between the input sequences accesible from the alignment property


    """
    aln = _PIR_alignment(msmseed.target_sequence, msmseed.target_id, msmseed.template_sequence, msmseed.template_id)
    msmseed.alignment = str(aln)
    return msmseed

def make_model(msmseed):
    """
    Use MODELLER from the Sali lab to create a model between the target and template specified in the input

    Parameters
    ----------
    msmseed : MSMSeed
        object containing the alignment between target and template and template structure

    Returns
    -------
    msmseed : MSMSeed
        object containing the homology model built from the input alignment and template structure
    """

    import tempfile
    import os
    import modeller
    import modeller.automodel
    import shutil
    import simtk.openmm.app as app
    #first, we need to make a temp directory where we can put the files MODELLER needs
    temp_dir = tempfile.mkdtemp()
    try:
        os.chdir(temp_dir)
        alignment_file = open('aln_tmp.pir','w')
        alignment_file.writelines(msmseed.alignment)
        alignment_file.close()
        template_file = open(msmseed.template_id + '.pdb','w')
        template_pdb = msmseed.template_structure
        template_pdb.writeFile(template_pdb.topology, template_pdb.positions, template_file)
        template_file.close()
        modeller.log.none()
        env = modeller.environ()
        env.io.atom_files_directory = temp_dir
        a = modeller.automodel.allhmodel(env,
                                         # file with template codes and target sequence
                                         alnfile  = 'aln_tmp.pir',
                                         # PDB codes of the template
                                         knowns   = msmseed.template_id,
                                         # code of the target
                                         sequence = msmseed.target_id)
        a.make()
        tmp_model_pdbfilename = a.outputs[0]['name']
        msmseed.target_model = app.PDBFile(tmp_model_pdbfilename)
        msmseed.target_restraints = open('%s.rsr' % msmseed.target_id, 'r').readlines()
    except:
        msmseed.error_message = 'MSMSeeder failed at the modelling stage'
        msmseed.error_state = -2
    finally:
        shutil.rmtree(temp_dir)
    return msmseed

def refine_implicitMD(msmseed, openmm_platform='CPU', niterations=100, nsteps_per_iteration=5):
    """
    Use OpenMM to perform an implicit solvent calculation on the output of MODELLER.

    Parameters
    ----------
    msmseed : MSMSeed
        object containing the modeled structure to be used for implicit simulation
    openmm_platofm : String, optional
        name of the OpenMM platform to be used in the simulation. Default CPU.
    niterations : int, optional
        number of iterations of steps to perform when running dynamics. Default 100.
    nsteps_per_iteration : int, optional
        number of steps to take each iteration. Default 5



    """

    import simtk.openmm as openmm
    import simtk.unit as unit
    import simtk.openmm.app as app

    forcefields_to_use = ['amber99sbildn.xml', 'amber99_obc.xml'] # list of forcefields to use in parameterization

    timestep = 2.0 * unit.femtoseconds # timestep
    temperature = 300.0 * unit.kelvin # simulation temperature
    collision_rate = 20.0 / unit.picoseconds # Langevin collision rate
    cutoff = None # nonbonded cutoff
    minimization_tolerance = 10.0 * unit.kilojoules_per_mole / unit.nanometer
    minimization_steps = 20

    kB = unit.MOLAR_GAS_CONSTANT_R
    kT = kB * temperature

    pH = 8.0

    forcefield = app.ForceField(*forcefields_to_use)
    platform = openmm.Platform.getPlatformByName(openmm_platform)
        #this will just be None if there is no gpu
    gpuid = os.getenv("CUDA_VISIBLE_DEVICES")
    if openmm_platform == 'CUDA':
        #here, gpuid is returned as a string, so it can be directly given to the platform method
        platform.setPropertyDefaultValue('CudaDeviceIndex', gpuid)
    if openmm_platform == 'OpenCL':
        platform.setPropertyDefaultValue('OpenCLDeviceIndex', gpuid)
    modeller = app.Modeller(msmseed.target_model.topology, msmseed.target_model.positions)
    modeller.addHydrogens(forcefield, pH=pH, variants=None)
    topology = modeller.getTopology()
    positions = modeller.getPositions()

    system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)

    openmm.LocalEnergyMinimizer.minimize(context, minimization_tolerance, minimization_steps)

    for iteration in range(niterations):
            # integrate dynamics
            integrator.step(nsteps_per_iteration)
            # get current state
            state = context.getState(getEnergy=True, getPositions=True)
            simulation_time = state.getTime()
            potential_energy = state.getPotentialEnergy()
            kinetic_energy = state.getKineticEnergy()
            import numpy
            if numpy.isnan(potential_energy/kT) or numpy.isnan(kinetic_energy/kT):
                msmseed.error_message = "Potential or kinetic energies are nan."
                msmseed.error_state=-2
                return msmseed
    state = context.getState(getPositions=True)
    refined = MDSys(modeller.topology, state.getPositions())
    msmseed.implicit_refined_model = refined
    return msmseed



def solvate_models(msmseed):
    """
    Calculate the number of waters needed to solvate the implicitly-refined model

    Parameters
    ----------
    msmseed : MSMSeed
        object containing the unsolvated, implicitly refined model structure

    Returns
    -------
    msmseed : MSMSeed
        object additionally containing the number of waters required to solvate the model

    """
    import os
    import simtk.unit as unit
    import simtk.openmm.app as app


    forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml'] # list of forcefields to use in parameterization
    nparticles_per_water = 3
    padding = 10.0 * unit.angstroms

    forcefield = app.ForceField(*forcefields_to_use)
    topology = msmseed.implicit_refined_model.topology
    positions = msmseed.implicit_refined_model.positions
    natoms_initial = len(positions)
    modeller = app.Modeller(topology, positions)
    modeller.addSolvent(forcefield, model='tip3p', padding=padding)
    positions = modeller.getPositions()
    natoms_final = len(positions)
    nwaters = (natoms_final - natoms_initial) / nparticles_per_water
    msmseed.nwaters = nwaters
    return msmseed



def calculate_nwaters(nwaters_list):
    """
    Calculate the number of waters at the 68th percentile.

    Parameters
    ----------
    nwaters_list : list of ints
        list of the number of waters needed to solvate each model

    Returns
    -------
    index68 : int
        the target number of waters at the 68th percentile
    """
    import numpy as np
    nwaters_array = np.array(nwaters_list)
    nwaters_array.sort()
    index68 = int((len(nwaters_array) - 1) * 0.68)
    #index95 = int((len(nwaters_array) - 1) * 0.95)
    return nwaters_array[index68]

def solvate_models_to_target(msmseed, target_nwaters):
    """
    Solvate the model to the target number of waters. If the model requires more than target_nwaters for solvation,
    it will be rejected.

    Parameters
    ----------
    msmseed : MSMSeed
        object containing an implicitly refined model to solvate
    target_nwaters : int
        number of waters to add

    Returns
    -------
    msmseed : MSMSeed
        object containing a solvated model with target_nwaters, or containing an error message

    """

    import simtk.openmm.app as app
    import simtk.unit as unit
    natoms_per_solvent = 3
    refined_model = msmseed.implicit_refined_model
    # Count initial atoms.
    natoms_initial = len(refined_model.positions)
    forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml']
    model='tip3p'
    forcefield = app.ForceField(*forcefields_to_use)

    # Solvate with zero padding to determine min number of waters and minimal unit cell dimensions.
    modeller = app.Modeller(refined_model.topology, refined_model.positions)
    modeller.addSolvent(forcefield, model=model, padding=0.0*unit.angstroms)
    topology = modeller.getTopology()
    positions = modeller.getPositions()
    box_min = topology.getUnitCellDimensions()
    natoms_min = len(positions) # minimal number of atoms
    nwaters_min = (natoms_min - natoms_initial) / natoms_per_solvent # minimal number of waters
    volume_min = box_min[0] * box_min[1] * box_min[2]
    residues = [ r for r in topology.residues() ] # build a list of residues
    nresidues_min = len(residues) # number of residues
    if nwaters_min > target_nwaters:
        msmseed.error_state = -4
        msmseed.error_message = "The minimally solvated model has more than the target number of waters"
        return msmseed

    #Make a slightly enlarged box
    scale = 1.1
    modeller = app.Modeller(refined_model.topology, refined_model.positions)
    topology = modeller.getTopology()
    topology.setUnitCellDimensions(box_min * scale)
    modeller.addSolvent(forcefield, model=model)
    positions = modeller.getPositions()
    box_enlarged = topology.getUnitCellDimensions()
    natoms_enlarged = len(positions) # minimal number of atoms
    nwaters_enlarged = (natoms_enlarged - natoms_initial) / natoms_per_solvent # minimal number of waters
    volume_enlarged = box_enlarged[0] * box_enlarged[1] * box_enlarged[2]
    density = (nwaters_enlarged - nwaters_min) / (volume_enlarged - volume_min)
    over_target = False
    extra_nwaters = 100
    while not over_target:
        delta_volume = (target_nwaters + extra_nwaters - nwaters_min) / density
        scale = ((volume_min + delta_volume) / volume_min)**(1.0/3.0)
        delta_volume = (target_nwaters + extra_nwaters - nwaters_min) / density
        modeller = app.Modeller(refined_model.topology, refined_model.positions)
        topology = modeller.getTopology()
        topology.setUnitCellDimensions(box_min * scale)
        modeller.addSolvent(forcefield, model=model)
        positions = modeller.getPositions()
        topology = modeller.getTopology()
        natoms = len(positions) # minimal number of atoms
        nwaters = (natoms - natoms_initial) / natoms_per_solvent # minimal number of waters
        if (nwaters > target_nwaters):
            over_target = True
        else:
            extra_nwaters += 100

    # Delete waters to achieve target.
    ndelete = nwaters - target_nwaters
    if (ndelete > 0):
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
    solvated_mdsys = MDSys(topology, positions)
    msmseed.solvated_model = solvated_mdsys
    return msmseed









def refine_explicitMD(msmseed, openmm_platform='CPU', niterations=1, nsteps_per_iteration=5):
    """
    Run an explicit-solvent MD refinement on a solvated model structure.

    Parameters
    ----------
    msmseed : MSMSeed
        object containing solvated model to simulate
    openmm_platform : String, optional
        name of the OpenMM platform to use in the simulation. Default CPU
    niterations : int, optional
        number of iterations of integrator steps to take. Default 1.
    nsteps_per_iteration : int, optional
        number of steps to take per iteration. Default 5.

    Returns
    -------
    msmseed : MSMSeed
         object containing gzipped explicitly-refined model pdb, integrator xml, system xml, and state xml

    """
    import simtk.openmm as openmm
    import simtk.openmm.app as app
    import simtk.unit as unit
    import time
    import StringIO
    import gzip

    platform = openmm.Platform.getPlatformByName(openmm_platform)
    #this will just be None if there is no gpu
    gpuid = os.getenv("CUDA_VISIBLE_DEVICES")
    if openmm_platform == 'CUDA':
        #here, gpuid is returned as a string, so it can be directly given to the platform method
        platform.setPropertyDefaultValue('CudaDeviceIndex', gpuid)
    if openmm_platform == 'OpenCL':
        platform.setPropertyDefaultValue('OpenCLDeviceIndex', gpuid)

    forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml'] # list of forcefields to use in parameterization

    timestep = 2.0 * unit.femtoseconds # timestep
    temperature = 300.0 * unit.kelvin # simulation temperature
    pressure = 1.0 * unit.atmospheres # simulation pressure
    collision_rate = 20.0 / unit.picoseconds # Langevin collision rate
    barostat_period = 50
    niterations = 100 # number of iterations

    nonbondedMethod = app.PME

    minimization_tolerance = 10.0 * unit.kilojoules_per_mole / unit.nanometer
    minimization_steps = 20

    kB = unit.MOLAR_GAS_CONSTANT_R
    kT = kB * temperature

    forcefield = app.ForceField(*forcefields_to_use)
    solvated_model = msmseed.solvated_model
    system = forcefield.createSystem(solvated_model.topology, nonbondedMethod=nonbondedMethod, constraints=app.HBonds)
    barostat = openmm.MonteCarloBarostat(pressure, temperature, barostat_period)
    system.addForce(barostat)
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(solvated_model.positions)
    openmm.LocalEnergyMinimizer.minimize(context, minimization_tolerance, minimization_steps)
    context.setVelocitiesToTemperature(temperature)
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
    state = context.getState(getPositions=True)

    #save the pdb of the current model
    explicit_model_pdb = StringIO.StringIO()
    with gzip.GzipFile(fileobj = explicit_model_pdb, mode = 'w') as output:
        app.PDBFile.writeHeader(solvated_model.topology, file=output)
        app.PDBFile.writeModel(solvated_model.topology, state.getPositions(), file=output)
        app.PDBFile.writeFooter(solvated_model.topology, file=output)
    msmseed.explicit_refined_pdb = explicit_model_pdb.getvalue()

    #save the state
    serialized_state_stringio = StringIO.StringIO()
    state_xml = openmm.XmlSerializer.serialize(state)
    with gzip.GzipFile(fileobj=serialized_state_stringio, mode = 'w') as output:
         output.write(state_xml)
    msmseed.explicit_refined_state = serialized_state_stringio.getvalue()

    #save the integrator
    serialized_integrator_stringio = StringIO.StringIO()
    integrator_xml = openmm.XmlSerializer.serialize(integrator)
    with gzip.GzipFile(fileobj = serialized_integrator_stringio, mode = 'w') as buffer:
        buffer.write(integrator_xml)
    msmseed.explicit_refined_integrator = serialized_state_stringio.getvalue()

    #save the system
    serialized_system_stringio = StringIO.StringIO()
    system_xml = openmm.XmlSerializer.serialize(system)
    with gzip.GzipFile(fileobj = serialized_system_stringio, mode = 'w') as buffer:
       buffer.write(system_xml)
    msmseed.explicit_refined_system = serialized_state_stringio.getvalue()


    return msmseed










if __name__=="__main__":
    import os
    import simtk.openmm.app as app
    os.chdir("/Users/grinawap/musashi_modelling")
    #os.environ['PYTHONPATH']='/Library/modeller-9.13/modlib:'

    #get ready to model
    target_sequence = "".join(open('rrm2.fasta','r').readlines())
    template_sequence = "".join(open('rrm1.fasta','r').readlines())
    template_structure = app.PDBFile("RRM1.pdb")
    test_msm_seed = distributed.MSMSeed(target_sequence, template_sequence, template_structure)
    print test_msm_seed.template_sequence
    print test_msm_seed.target_sequence
    make_PIR_alignment(test_msm_seed)
    print test_msm_seed.alignment
    make_model(test_msm_seed)
    refine_implicitMD(test_msm_seed,niterations=1)
