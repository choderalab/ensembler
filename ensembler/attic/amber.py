# Module to facilitate simulating an OpenMM system in AMBER.

#
# IMPORTS
#

def get_last_good_chunk(model_dir, destructive=False, verbose=False):
    '''
    For a given model dir, looks for the most recent chunk which is in good condition (currently this just checks quality of the .rst file).
    Returns the chunk index.
    If no good chunks were found, the returned chunk index is -1.
    If destructive option is set, then any bad chunks found during this process (e.g. with truncated .rst file) are deleted.

    Requires amber.crd to be present in the model_dir - necessary to get the natoms, which is used to validate .rst file length.
    '''
    import os, shutil
    # First make a list of all previous chunks
    dirs_in_model_dir = os.walk(model_dir).next()[1]
    prev_chunk_indices = []
    for dir in dirs_in_model_dir:
        if dir[0:4] == 'traj':
            prev_chunk_indices.append(int(dir[4:]))

    last_good_chunk_index = -1

    if verbose: print('prev_chunk_indices for model_dir %s:' % model_dir, prev_chunk_indices)

    # If previous chunks exist, iterate through them in reverse order until finding one with a .rst file in good condition
    if len(prev_chunk_indices) > 0:
        prev_chunk_indices.sort(reverse=True)
        model_natoms = get_natoms(model_dir)

        for c in prev_chunk_indices:
            last_chunk_dir = os.path.join(model_dir, 'traj' + str(c))
            last_rst_file_good = test_rst_file(last_chunk_dir, model_natoms, verbose=verbose)
            if last_rst_file_good:
                last_good_chunk_index = c
                break
            # If the .rst file is not in good condition, and if the destructive option is set, delete the associated chunk
            else:
                if destructive == True:
                    shutil.rmtree(last_chunk_dir)
    # If a chunk with a good .rst file has not been found, then each chunk will subsequently have been deleted, and last_chunk_index will remain at -1

    return last_good_chunk_index


def write_rst_from_nc(traj_filepath='md.nc', parm_filepath='../amber.prmtop', rst_filepath='../md.rst', cpptraj_input_filepath='rst_from_nc.in', frame='last'):
    '''
    Caution: output .rst file will not contain velocities if the input trajectory does not.
    '''
    import subprocess
    # Create input file for cpptraj
    with open(cpptraj_input_filepath, 'w') as cpptraj_input_file:
        cpptraj_input_file.write('trajin %s lastframe\n' % traj_filepath)
        cpptraj_input_file.write('trajout %s\n' % rst_filepath)

    # Run cpptraj
    # NOTE may need to specify AMBERHOME - get this from sys.path
    subprocess.check_output(['cpptraj', '-i', cpptraj_input_filepath, '-p', parm_filepath])

def get_simlength(chunk_dir):
    '''
    For a given traj chunk, returns the cumulative traj length from the .rst file.
    '''
    import os
    rst_file_path = os.path.join(chunk_dir, 'md.rst')
    with open(rst_file_path, 'r') as rst_file:
        rst_file.readline()
        simlength = float(rst_file.readline().split()[-1].strip())
    return simlength

def get_natoms(model_dir):
    '''
    For a given model, gets the number of atoms from the input .crd file.
    '''
    import os
    input_coord_file_path = os.path.join(model_dir, 'amber.crd')
    with open(input_coord_file_path, 'r') as input_coord_file:
        input_coord_file.readline()
        natoms = int(input_coord_file.readline().strip())
    return natoms

def test_rst_file(chunk_dir, natoms, verbose=False):
    '''
    For a given model with a given natoms, carries out the following checks on the .rst file:
      * checks the .rst file exists
      * checks the number of lines in the .rst file is as expected from the natoms.
    Returns True if the .rst file is ok, False if not (and prints a wanring if verbose is set.)
    '''
    import os, subprocess
    rst_file_path = os.path.join(chunk_dir, 'md.rst')

    if not os.path.exists(rst_file_path):
        if verbose:
            print('WARNING for chunk_dir %s: file "md.rst" not found.' % chunk_dir)
        return False

    rst_nlines = int( subprocess.check_output(['wc', '-l', rst_file_path]).split()[0] )

    expected_rst_nlines = (((natoms / 2) + (natoms % 2)) * 2) + 3

    if rst_nlines != expected_rst_nlines:
        if verbose:
            print('WARNING for chunk_dir %s: only %d lines found in md.rst - expected %d lines from the number of atoms listed in amber.crd' % (chunk_dir, rst_nlines, expected_rst_nlines))
        return False
    else:
        return True

def _readFileContents(filename):
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
    
def _writeFileContents(filename, contents):
    outfile = open(filename, 'w')
    outfile.write(contents)
    outfile.close()

def createAmberInputFiles(topology, system, state, nproteinatoms, openmm_forcefields_to_use, amber_forcefield_to_use, prmtop_filename, inpcrd_filename, verbose=False, shell='/bin/tcsh'):
    '''
    Create AMBER input files from an OpenMM system.

    ARGUMENTS

    topology (simtk.openmm.app.Topology) - the topology of the system to be simulated
    system (simtk.openmm.System) - the System object for the system to be simulated
    state (simtk.openmm.State) - the State of the system to be simulated
    nproteinatoms (int) - number of protein atoms (required for centering)
    openmm_forcefields_to_use (list of string) - list of forcefield XML files to use in generating AMBER translation names
    amber_forcefield_to_use (string) - protein force field to be used in LEaP
    prmtop_filename (string) - the filename of the AMBER prmtop file to be created
    inpcrd_filename (string) - the filename of the AMBER inpcrd file to be created

    OPTIONAL ARGUMENTS

    verbose (boolean) - if True, will print verbose output

    NOTES

    Note that AmberTools must be installed.  LEaP is used to create the input files.

    TODO

    * ?Use a single forcefields_to_use argument to automatically select the appropriate AMBER and OpenMM forcefield files.
    * Use temporary directory and copy out to desired filenames.
    * Decide on whether to keep ambernames.pdb or amber.pdb in each directory

    RETURNS
    
    None

    '''    

    from simtk import openmm
    from simtk import unit
    
    # Create box and center protein
    if verbose: print("Creating box and centering protein in unit cell...")
    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName('Reference')
    context = openmm.Context(system, integrator, platform)
    context.setState(state)
    state = context.getState(getPositions=True, enforcePeriodicBox=True)
    positions = state.getPositions(asNumpy=True)
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
    mean = unit.Quantity((positions[0:nproteinatoms,:] / unit.angstrom).mean(0), unit.angstroms)
    for i in range(system.getNumParticles()):
        positions[i,:] -= mean[:] + box_vectors[0,0]/2.0
    context.setPositions(positions)
    state = context.getState(getPositions=True, enforcePeriodicBox=True)
    positions = state.getPositions(asNumpy=True)
    del context

    # Configure text for LEaP input file
    leap_template = '''
# Set up solvated protein system for explicit solvent simulation.

# Load AMBER ff99sb-ildn forcefield for protein.
source %(amber_forcefield_to_use)s

# Load modified ion parameters.
loadAmberParams frcmod.ionsjc_tip3p

# Load in protein (with all residues modeled in).
system = loadPdb ambernames.pdb

# Generate box.
solvatebox system TIP3PBOX 0.0001 10000 iso

# Check protein.
check system

# Report on net charge.
charge system

# Write parameters.
saveAmberParm system amber.prmtop amber.crd

# Exit
quit
''' % vars()

    # Write LEaP input file.
    if verbose: print("Writing LEaP input file...")
    leap_filename = 'setup.leap.in'
    outfile = open(leap_filename, 'w')
    outfile.write(leap_template)
    outfile.close()

    # Clear leap.log.
    import os, os.path
    leap_log_filename = 'leap.log'
    if os.path.exists(leap_log_filename):
        os.remove(leap_log_filename)
    
    # Determine atom names.
    (atoms, sorted_atom_indices) = _assignNamesFromForceFieldTemplates(topology, system, openmm_forcefields_to_use)
    
    # Re-sort atoms.
    resort_atoms = True
    if resort_atoms:
        import copy
        atoms2 = copy.deepcopy(atoms)
        positions2 = copy.deepcopy(positions)
        for (new_index, old_index) in enumerate(sorted_atom_indices):
            atoms2[new_index] = atoms[old_index]
            atoms2[new_index].index = new_index+1
            positions2[new_index,:] = positions[old_index,:]
    
        atoms = atoms2
        positions = positions2

    # Write PDB file with AMBER names for atoms and residues.
    amberpdb_filename = 'ambernames.pdb'
    outfile = open(amberpdb_filename, 'w')
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
    index = 0
    wrap_coordinates = False
    for atom in atoms:
        residue = atom.residue
        atomIndex = atom.index
        atomName = atom.name
        resName = atom.resname
        chainName = chr(ord('A') + residue.chain.index)
        resIndex = residue.index

        if len(resName) > 3:
            resName = resName[1:]
        outfile.write("ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n" % (atomIndex%100000, atomName, resName, chainName, (resIndex+1)%10000, positions[index,0]/unit.angstrom, positions[index,1]/unit.angstrom, positions[index,2]/unit.angstrom))
        index += 1
    outfile.close()

    # Run tleap.
    import subprocess
    command = 'setenv AMBERHOME $AMBERHOME_CPU; tleap -f setup.leap.in >& setup.leap.out'
    try:
        output = subprocess.check_output(command, shell=True, executable=shell, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print(e, 'with output:\n' + e.output)
        raise e
    if verbose: print(output)

    # Write periodic AMBER .crd file.
    overwrite_crd = True
    if overwrite_crd:
        inpcrd_filename = 'amber.crd'
        outfile = open(inpcrd_filename, 'w')
        outfile.write('Automatically converted from OpenMM\n')
        natoms = system.getNumParticles()
        outfile.write('%6d\n' % natoms)
        nwritten = 0
        for i in range(natoms):
            if nwritten == 6:
                outfile.write('\n')
                nwritten = 0
            for k in range(3):
                outfile.write('%12.7f' % (positions[i][k] / unit.angstroms))
                nwritten += 1
        outfile.write('\n')
        for k in range(3):
            outfile.write('%12.7f' % (box_vectors[k][k]/unit.angstroms))
        for k in range(3):
            outfile.write('%12.7f' % 90.0)
        outfile.write('\n')
        outfile.close()

        # Generate PDB file.
        import subprocess
        command = 'cat amber.crd | ambpdb -p amber.prmtop -aatm > amber.pdb'
        try:
            output = subprocess.check_output(command, shell=True, executable=shell, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(e, 'with output:\n' + e.output)
            raise e
        print(output)

    return

def _assignNamesFromForceFieldTemplates(topology, system, openmm_forcefields_to_use):
    '''
    Assign Amber atom and residue names from specified forcefield templates.
    Requires OpenMM version >= 5.2.

    ARGUMENTS

    topology
    system
    openmm_forcefields_to_use (list of strings) - list of XML files of forcefields containing templates

    RETURNS
    
    atoms - list of atom information 
    sorted_atom_indices (list of int) - atom indices in order that they appear in forcefield templates, in case re-sorting is required

    '''

    # Assign AMBER atom and residue names.
    from simtk.openmm.app import ForceField
    import simtk.openmm.app.forcefield
    forcefield = ForceField(*openmm_forcefields_to_use)
        
    data = ForceField._SystemData()
    atomIndices = {}
    for index, atom in enumerate(topology.atoms()):
        data.atoms.append(atom)
        atomIndices[atom] = index

    # Make a list of all bonds
    for bond in topology.bonds():
        if bond[0] in atomIndices and bond[1] in atomIndices:
            data.bonds.append(ForceField._BondData(atomIndices[bond[0]], atomIndices[bond[1]]))

    # Record which atoms are bonded to each other atom
    bondedToAtom = []
    for i in range(len(data.atoms)):
        bondedToAtom.append(set())
        data.atomBonds.append([])
    for i in range(len(data.bonds)):
        bond = data.bonds[i]
        bondedToAtom[bond.atom1].add(bond.atom2)
        bondedToAtom[bond.atom2].add(bond.atom1)
        data.atomBonds[bond.atom1].append(i)
        data.atomBonds[bond.atom2].append(i)
        
    # Find the template matching each residue and assign atom types.        
    sorted_atom_indices = list()
    for chain in topology.chains():
        for res in chain.residues():
            template = None
            matches = None
            # NOTE: Before OpenMM 5.2, was necessary to convert output of _createResidueSignature using the method simtk.openmm.app.forcefield._signatureToString - this is now done within _createResidueSignature
            sig = simtk.openmm.app.forcefield._createResidueSignature([atom.element for atom in res.atoms()])
            if sig != '':
                if sig in forcefield._templateSignatures:
                    for t in forcefield._templateSignatures[sig]:
                        # NOTE: No longer necessary to pass atomIndices in OpenMM > 5.2
                        matches = simtk.openmm.app.forcefield._matchResidue(res, t, bondedToAtom)
                        if matches is not None:
                            template = t
                            break
            if matches is None:
                # Check templates involving virtual sites
                for t in forcefield._templateSignatures[None]:
                    matches = simtk.openmm.app.forcefield._matchResidue(res, t, bondedToAtom)
                    if matches is not None:
                        template = t
                        break
            if matches is None:
                raise ValueError('No template found for residue %d (%s).  This might mean your input topology is missing some atoms or bonds, or possibly that you are using the wrong force field.' % (res.index+1, res.name))

            # Sort matches by order in template.
            atom_indices = [ atom.index for atom in res.atoms() ]
            for local_index in range(len(template.atoms)):
                if local_index in matches:
                    sorted_atom_indices.append(min(atom_indices) + matches.index(local_index))

            for (atom, match) in zip(res.atoms(), matches):
                data.atomType[atom] = template.atoms[match].type
                # Rename atom (JDC).
                atom.name = template.atoms[match].name
                atom.resname = template.name

                for site in template.virtualSites:
                    if match == site.index:
                        data.virtualSites[atom] = site

    return (data.atoms, sorted_atom_indices)

def create_mdin(ntimesteps, restart=False, coord_writeperiod=50000, energy_writeperiod=50000, rst_writeperiod=50000, force=False):
    '''
    restart option creates an mdin file appropriate for restarting a simulation
    force option will overwrite an existing mdin file

    TODO:

    Allow all parameters to be specified with OpenMM-like interface?
    
    '''
    import os
    # Return if file already exists (and if not forcing overwrite)
    if restart:
        mdin_filename = 'md.restart.in'
    else:
        mdin_filename = 'md.start.in'
    if os.path.exists(mdin_filename) and force == False:
        return

    # Set up some parameters based on whether or not this simulation is being restarted
    if restart:
        irest=1
        ntx=5
    else:
        irest=0
        ntx=0

    mdin_template = '''\
Langevin dynamics in periodic box.
 &cntrl
   imin=0, irest=%(irest)d, ntx=%(ntx)d,
   imin=0, irest=0, ntx=1,
   ntc=2, ntf=2, 
   nstlim=%(ntimesteps)d, 
   ntpr=%(energy_writeperiod)d, ntwx=%(coord_writeperiod)d,
   ntwr=%(rst_writeperiod)d, 
   dt=0.002, cut=9.,
   ntt=3, gamma_ln=5.0, tempi=300.0, temp0=300.0,
   ntb=2, ntp=1, pres0=1.01317122594, taup=10.0,
   ioutfm=1,
 &end
''' % vars()
    outfile = open(mdin_filename, 'w')
    outfile.write(mdin_template)
    outfile.close()
                        
def run_sander(verbose=False, restart=False, overwrite_mdin=False):
    # Generate MD input file.
    create_mdin(ntimesteps, restart=restart)

    # Run sander.
    import deprecated_commands
    print("Running sander...")
    command = 'sander -O -i md.sander.in -o md.sander.out -p amber.prmtop -c amber.crd -r md.rst -x md.nc -e md.ene'
    output = deprecated_commands.getoutput(command)
    print(output)

    return

def run_pmemd(ncpus=1, shell='/bin/tcsh', verbose=False, restart=False, overwrite_mdin=False):
    # Generate MD input file.
    create_mdin(ntimesteps, restart=restart)

    # Run CPU version of PMEMD
    import subprocess
    command = "setenv AMBERHOME $AMBERHOME_CPU; $AMBERHOME/bin/pmemd -O -i md.sander.in -o md.sander.out -p amber.prmtop -c amber.crd -inf amber.mdinfo -x amber.nc -r amber.restrt" % vars() 
    output = subprocess.check_output(command, shell=True, executable=shell)
    if verbose: print(output)

    return

def run_pmemd_cuda(ntimesteps, gpuid=1, shell='/bin/tcsh', verbose=False, restart=False, overwrite_mdin=False):
    '''
Run serial GPU version of PMEMD.
One CPU is used per GPU.

This routine is designed to run in a chunk dir (named "traj0", "traj1", ...).
Topology and coord input are taken from the parent model dir.
MD simulation params are taken from the chunk dir.
Trajectory, simulation info and restart file are output to each the dir.

If restarting a trajectory, set restart to the name of the restart file to be used as input.
    '''
    import subprocess
    cuda_visible_devices = gpuid

    # Generate MD input file
    create_mdin(ntimesteps, restart=restart)

    # Run Serial CUDA version of PMEMD.
    if restart:    
        command = "setenv AMBERHOME $AMBERHOME_GPU; setenv CUDA_VISIBLE_DEVICES %(cuda_visible_devices)s; $AMBERHOME/bin/pmemd.cuda -O -i md.restart.in -p ../amber.prmtop -c %(restart)s -o md.out -inf md.mdinfo -x md.nc -r md.rst" % vars() 
    else:
        command = "setenv AMBERHOME $AMBERHOME_GPU; setenv CUDA_VISIBLE_DEVICES %(cuda_visible_devices)s; $AMBERHOME/bin/pmemd.cuda -O -i md.start.in -p ../amber.prmtop -c ../amber.crd -o md.out -inf md.mdinfo -x md.nc -r md.rst" % vars() 

    output = subprocess.check_output(command, shell=True, executable=shell)
    if verbose: print(output)

    return

def run_pmemd_cuda_mpi(gpupn=1, shell='/bin/tcsh', verbose=False, restart=False, overwrite_mdin=False):
    '''
Run parallel GPU version of PMEMD.
One CPU used per node.
TODO
    '''
    import subprocess
    cuda_visible_devices = ''
    for i in range(gpupn):
	if (i > 0): cuda_visible_devices += ','
	cuda_visible_devices += '%d' % i

    return

#
# TEST HARNESS
#

if __name__ == '__main__':

    # Read OpenMM System and coordinates.
    from simtk import openmm
    from simtk import unit
    system = openmm.XmlSerializer.deserialize(_readFileContents('explicit-system.xml'))
    state = openmm.XmlSerializer.deserialize(_readFileContents('explicit-state.xml'))

    # Read explicitly solvated PDB file.
    from simtk.openmm.app import PDBFile
    pdb_filename = 'explicit-refined.pdb'
    pdbfile = PDBFile(pdb_filename)
    topology = pdbfile.getTopology()

    # Create AMBER input files.
    prmtop_filename = 'amber.prmtop'
    inpcrd_filename = 'amber.crd'
    openmm_forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml']
    amber_forcefield_to_use = '/mnt/b/projects/sciteam/jn6/GIT/amber-gnu/dat/leap/cmd/oldff/leaprc.ff99SBildn'
    createAmberInputFiles(topology, system, state, openmm_forcefields_to_use, amber_forcefield_to_use, prmtop_filename, inpcrd_filename, verbose=True)

    # Run dynamics.
    ntimesteps = 5000
    run_pmemd_cuda(ntimesteps, verbose=True)

