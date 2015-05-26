import os
import datetime
import traceback
import gzip
import sys
import subprocess
import yaml
from yaml.scanner import ScannerError
import warnings
import socket
from collections import deque
import numpy as np
import Bio
import ensembler
import ensembler.version
from ensembler.core import mpistate, logger
import simtk.unit as unit
import simtk.openmm as openmm
import simtk.openmm.app as app
import simtk.openmm.version


def refine_implicit_md(
        openmm_platform=None, gpupn=1, process_only_these_targets=None,
        process_only_these_templates=None, template_seqid_cutoff=None,
        verbose=False, write_trajectory=False,
        include_disulfide_bonds=False,
        ff='amber99sbildn',
        implicit_water_model='amber99_obc',
        sim_length=100.0 * unit.picoseconds,
        timestep=2.0 * unit.femtoseconds,             # timestep
        temperature=300.0 * unit.kelvin,              # simulation temperature
        collision_rate=20.0 / unit.picoseconds,       # Langevin collision rate
        cutoff=None,                                  # nonbonded cutoff
        minimization_tolerance=10.0 * unit.kilojoules_per_mole / unit.nanometer,
        minimization_steps=20,
        nsteps_per_iteration=500,
        ph=7.0,
        retry_failed_runs=False,
        cpu_platform_threads=1):
    # TODO - refactor
    '''Run MD refinement in implicit solvent.

    MPI-enabled.
    '''
    gpuid = mpistate.rank % gpupn

    models_dir = os.path.abspath(ensembler.core.default_project_dirnames.models)

    targets, templates_resolved_seq = ensembler.core.get_targets_and_templates()

    if process_only_these_templates:
        selected_template_indices = [i for i, seq in enumerate(templates_resolved_seq) if seq.id in process_only_these_templates]
    else:
        selected_template_indices = range(len(templates_resolved_seq))

    if not openmm_platform:
        openmm_platform = auto_select_openmm_platform()

    if openmm_platform == 'CPU':
        platform_properties = {'CpuThreads': str(cpu_platform_threads)}
    else:
        platform_properties = {}

    ff_files = [ff+'.xml', implicit_water_model+'.xml']
    forcefield = app.ForceField(*ff_files)

    kB = unit.MOLAR_GAS_CONSTANT_R
    kT = kB * temperature

    niterations = int((sim_length / timestep) / nsteps_per_iteration)

    def simulate_implicit_md():

        if verbose: print("Reading model...")
        with gzip.open(model_filename) as model_file:
            pdb = app.PDBFile(model_file)

        # Set up Platform
        platform = openmm.Platform.getPlatformByName(openmm_platform)
        if 'CUDA_VISIBLE_DEVICES' not in os.environ:
            # Set GPU id.
            if openmm_platform == 'CUDA':
                platform.setPropertyDefaultValue('CudaDeviceIndex', '%d' % gpuid)
            elif openmm_platform == 'OpenCL':
                platform.setPropertyDefaultValue('OpenCLDeviceIndex', '%d' % gpuid)

        # Construct Modeller object with same topology as ref structure
        # (necessary to keep disulfide bonds consistent)
        modeller = app.Modeller(reference_topology, pdb.positions)
        # set_openmm_topology_bonds_from_atom_indices(modeller.topology, reference_bonds)
        # Add missing protons.
        modeller.addHydrogens(forcefield, pH=ph, variants=reference_variants)
        topology = modeller.getTopology()
        positions = modeller.getPositions()

        if verbose: print("Constructing System object...")
        if cutoff is None:
            system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
        else:
            system = forcefield.createSystem(topology, nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=cutoff, constraints=app.HBonds)

        if verbose: print("Creating Context...")
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator, platform, platform_properties)
        context.setPositions(positions)

        if verbose: print("Minimizing structure...")
        openmm.LocalEnergyMinimizer.minimize(context, minimization_tolerance, minimization_steps)

        if write_trajectory:
            # Open trajectory for writing.
            if verbose: print("Opening trajectory for writing...")
            trajectory_filename = os.path.join(model_dir, 'implicit-trajectory.pdb.gz')
            trajectory_outfile = gzip.open(trajectory_filename, 'w')
            app.PDBFile.writeHeader(topology, file=trajectory_outfile)

        # Open energy trajectory for writing
        energy_filename = os.path.join(model_dir, 'implicit-energies.txt')
        energy_outfile = open(energy_filename, 'w')
        energy_outfile.write('# iteration | simulation time (ps) | potential_energy (kT) | kinetic_energy (kT) | ns per day\n')

        if verbose: print("Running dynamics...")
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
            if verbose: print(
                "  %8.1f ps : potential %8.3f kT | kinetic %8.3f kT | %.3f ns/day | %.3f s remain"
                % (
                    simulation_time / unit.picoseconds, potential_energy / kT, kinetic_energy / kT,
                    ns_per_day,
                    elapsed_time * (niterations-iteration-1) / (iteration+1) / unit.seconds
                )
            )

            # Check energies are still finite.
            if np.isnan(potential_energy/kT) or np.isnan(kinetic_energy/kT):
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
        pdb_outfile = gzip.open(pdb_filename, 'w')
        app.PDBFile.writeHeader(topology, file=pdb_outfile)
        app.PDBFile.writeFile(topology, state.getPositions(), file=pdb_outfile)
        app.PDBFile.writeFooter(topology, file=pdb_outfile)
        pdb_outfile.close()




    for target in targets:
        if process_only_these_targets and (target.id not in process_only_these_targets):
            continue
        models_target_dir = os.path.join(models_dir, target.id)
        if mpistate.rank == 0:
            target_starttime = datetime.datetime.utcnow()
            if not os.path.exists(models_target_dir):
                continue

        mpistate.comm.Barrier()

        # ========
        # Determine topology (including protonation state) to use throughout
        # ========

        reference_model_id = get_highest_seqid_existing_model(models_target_dir=models_target_dir)
        if reference_model_id is None:
            continue

        reference_model_path = os.path.join(models_target_dir, reference_model_id, 'model.pdb.gz')

        with gzip.open(reference_model_path) as reference_pdb_file:
            reference_pdb = app.PDBFile(reference_pdb_file)

        logger.debug("Using %s as highest identity model" % (reference_model_id))

        if not include_disulfide_bonds:
            remove_disulfide_bonds_from_topology(reference_pdb.topology)

        # Build topology for reference model
        modeller = app.Modeller(reference_pdb.topology, reference_pdb.positions)
        reference_topology = modeller.topology
        reference_variants = modeller.addHydrogens(forcefield, pH=ph)
        if verbose:
            print("Reference variants extracted:")
            if reference_variants != None:
                for (residue_index, residue) in enumerate(reference_variants):
                    if residue != None:
                        print("%8d %s" % (residue_index+1, residue))
                print("")
            else: print(reference_variants)

        if template_seqid_cutoff:
            process_only_these_templates = ensembler.core.select_templates_by_seqid_cutoff(target.id, seqid_cutoff=template_seqid_cutoff)
            selected_template_indices = [i for i, seq in enumerate(templates_resolved_seq) if seq.id in process_only_these_templates]

        ntemplates_selected = len(selected_template_indices)

        for template_index in range(mpistate.rank, ntemplates_selected, mpistate.size):
            template = templates_resolved_seq[selected_template_indices[template_index]]

            model_dir = os.path.join(models_target_dir, template.id)
            if not os.path.exists(model_dir): continue

            # Only simulate models that are unique following filtering by clustering.
            unique_by_clustering = os.path.exists(os.path.join(model_dir, 'unique_by_clustering'))
            if not unique_by_clustering: continue

            # Pass if this simulation has already been run.
            log_filepath = os.path.join(model_dir, 'implicit-log.yaml')
            if os.path.exists(log_filepath):
                with open(log_filepath) as log_file:
                    log_data = yaml.load(log_file, Loader=ensembler.core.YamlLoader)
                    if log_data.get('successful') is True:
                        continue
                    if log_data.get('finished') is True and (retry_failed_runs is False and log_data.get('successful') is False):
                        continue

            # Check to make sure the initial model file is present.
            model_filename = os.path.join(model_dir, 'model.pdb.gz')
            if not os.path.exists(model_filename):
                if verbose: print('model.pdb.gz not present: target %s template %s rank %d gpuid %d' % (target.id, template.id, mpistate.rank, gpuid))
                continue

            pdb_filename = os.path.join(model_dir, 'implicit-refined.pdb.gz')

            print("-------------------------------------------------------------------------")
            print("Simulating %s => %s in implicit solvent for %.1f ps (MPI rank: %d, GPU ID: %d)" % (target.id, template.id, niterations * nsteps_per_iteration * timestep / unit.picoseconds, mpistate.rank, gpuid))
            print("-------------------------------------------------------------------------")

            # Open log file
            log_data = {
                'mpi_rank': mpistate.rank,
                'gpuid': gpuid if 'CUDA_VISIBLE_DEVICES' not in os.environ else os.environ['CUDA_VISIBLE_DEVICES'],
                'openmm_platform': openmm_platform,
                'sim_length': '%s' % sim_length,
                'finished': False,
                }
            log_file = ensembler.core.LogFile(log_filepath)
            log_file.log(new_log_data=log_data)

            try:
                start = datetime.datetime.utcnow()
                simulate_implicit_md()
                timing = ensembler.core.strf_timedelta(datetime.datetime.utcnow() - start)
                log_data = {
                    'finished': True,
                    'timing': timing,
                    'successful': True,
                    }
                log_file.log(new_log_data=log_data)
            except Exception as e:
                trbk = traceback.format_exc()
                warnings.warn(
                    '= ERROR start: MPI rank {0} hostname {1} gpuid {2} =\n{3}\n{4}\n= ERROR end: MPI rank {0} hostname {1} gpuid {2}'.format(
                        mpistate.rank, socket.gethostname(), gpuid, e, trbk
                    )
                )
                timing = ensembler.core.strf_timedelta(datetime.datetime.utcnow() - start)
                log_data = {
                    'exception': e,
                    'traceback': ensembler.core.literal_str(trbk),
                    'timing': timing,
                    'finished': True,
                    'successful': False,
                    }
                log_file.log(new_log_data=log_data)

        if verbose:
            print('Finished template loop: rank %d' % mpistate.rank)

        mpistate.comm.Barrier()

        if mpistate.rank == 0:
            project_metadata = ensembler.core.ProjectMetadata(project_stage='refine_implicit_md', target_id=target.id)

            datestamp = ensembler.core.get_utcnow_formatted()
            nsuccessful_refinements = subprocess.check_output(['find', models_target_dir, '-name', 'implicit-refined.pdb.gz']).count('\n')
            target_timedelta = datetime.datetime.utcnow() - target_starttime

            metadata = {
                'target_id': target.id,
                'datestamp': datestamp,
                'template_seqid_cutoff': template_seqid_cutoff,
                'process_only_these_targets': process_only_these_targets,
                'process_only_these_templates': process_only_these_templates,
                'timing': ensembler.core.strf_timedelta(target_timedelta),
                'ff': ff,
                'implicit_water_model': implicit_water_model,
                'nsuccessful_refinements': nsuccessful_refinements,
                'python_version': sys.version.split('|')[0].strip(),
                'python_full_version': ensembler.core.literal_str(sys.version),
                'ensembler_version': ensembler.version.short_version,
                'ensembler_commit': ensembler.version.git_revision,
                'biopython_version': Bio.__version__,
                'openmm_version': simtk.openmm.version.short_version,
                'openmm_commit': simtk.openmm.version.git_revision,
            }

            project_metadata.add_data(metadata)
            project_metadata.write()

        mpistate.comm.Barrier()

    mpistate.comm.Barrier()
    if mpistate.rank == 0:
        print('Done.')


def auto_select_openmm_platform():
    for platform_name in ['CUDA', 'OpenCL', 'CPU', 'Reference']:
        try:
            platform = openmm.Platform.getPlatformByName(platform_name)
            if type(platform) == openmm.Platform:
                logger.info('Auto-selected OpenMM platform: %s' % platform_name)
                return platform_name
        except Exception:
            continue
    raise Exception('No OpenMM platform found')


def get_highest_seqid_existing_model(targetid=None, models_target_dir=None):
    """
    Parameters
    ----------
    targetid: str
    models_target_dir: str

    Returns
    -------
    reference_model_id: str
        e.g. 'FAK1_HUMAN_D0_4KAB_B'
    """
    if not models_target_dir and targetid:
        models_target_dir = os.path.join(ensembler.core.default_project_dirnames.models, targetid)

    seqids_filepath = os.path.join(models_target_dir, 'sequence-identities.txt')
    if not os.path.exists(seqids_filepath):
        warnings.warn('ERROR: sequence-identities.txt file not found at path %s' % seqids_filepath)
        return None

    with open(seqids_filepath, 'r') as seqids_file:
        seqids_data = [line.split() for line in seqids_file.readlines()]

    # Find highest sequence identity model - topology will be used for all models
    for seqid_data in seqids_data:
        reference_model_id, reference_identity = seqid_data
        reference_pdb_filepath = os.path.join(models_target_dir, reference_model_id, 'model.pdb.gz')
        if os.path.exists(reference_pdb_filepath):
            return reference_model_id

    warnings.warn('ERROR: reference PDB model not found at path')
    return None


def remove_disulfide_bonds_from_topology(topology):
    """Should work with topology object from OpenMM or mdtraj.

    Parameters
    ----------
      topology: simtk.openmm.app.Topology or mdtraj.Topology
    """
    remove_bond_indices = []
    for b, bond in enumerate(topology._bonds):
        atom0, atom1 = bond
        if (
            atom0.residue.name == 'CYS' and atom1.residue.name == 'CYS'
            and (atom0.residue.index != atom1.residue.index)
            and (atom0.name == 'SG' and atom0.name == 'SG')
            ):
            remove_bond_indices.append(b)
    [topology._bonds.pop(b) for b in remove_bond_indices]


def solvate_models(process_only_these_targets=None, process_only_these_templates=None,
                   template_seqid_cutoff=None,
                   ff='amber99sbildn',
                   water_model='tip3p',
                   verbose=False,
                   padding=None):
    """Solvate models which have been subjected to MD refinement with implicit solvent.

    MPI-enabled.
    """
    if padding is None:
        padding = 10.0 * unit.angstroms
    elif type(padding) is float:
        padding = padding * unit.angstroms
    else:
        raise Exception('padding must be passed as a float (in Angstroms)')

    models_dir = os.path.abspath(ensembler.core.default_project_dirnames.models)

    targets, templates_resolved_seq = ensembler.core.get_targets_and_templates()

    if process_only_these_templates:
        selected_template_indices = [i for i, seq in enumerate(templates_resolved_seq) if seq.id in process_only_these_templates]
    else:
        selected_template_indices = range(len(templates_resolved_seq))

    ff_files = [ff+'.xml', water_model+'.xml']
    forcefield = app.ForceField(*ff_files)

    for target in targets:

        if process_only_these_targets and (target.id not in process_only_these_targets): continue

        models_target_dir = os.path.join(models_dir, target.id)
        if not os.path.exists(models_target_dir): continue

        if mpistate.rank == 0:
            target_starttime = datetime.datetime.utcnow()

        if template_seqid_cutoff:
            process_only_these_templates = ensembler.core.select_templates_by_seqid_cutoff(target.id, seqid_cutoff=template_seqid_cutoff)
            selected_template_indices = [i for i, seq in enumerate(templates_resolved_seq) if seq.id in process_only_these_templates]

        ntemplates_selected = len(selected_template_indices)

        for template_index in range(mpistate.rank, ntemplates_selected, mpistate.size):
            template = templates_resolved_seq[selected_template_indices[template_index]]

            model_dir = os.path.join(models_target_dir, template.id)
            if not os.path.exists(model_dir): continue

            model_filename = os.path.join(model_dir, 'implicit-refined.pdb.gz')
            if not os.path.exists(model_filename): continue

            print("-------------------------------------------------------------------------")
            print("Solvating %s => %s in explicit solvent" % (target.id, template.id))
            print("-------------------------------------------------------------------------")
            
            # Pass if solvation has already been run for this model.
            nwaters_filename = os.path.join(model_dir, 'nwaters.txt')
            if os.path.exists(nwaters_filename):
                continue

            try:
                if verbose: print("Reading model...")
                with gzip.open(model_filename) as model_file:
                    pdb = app.PDBFile(model_file)

                # Count initial atoms.
                natoms_initial = len(pdb.positions)

                # Add solvent
                if verbose: print("Solvating model...")
                modeller = app.Modeller(pdb.topology, pdb.positions)
                modeller.addSolvent(forcefield, model='tip3p', padding=padding)
                positions = modeller.getPositions()

                # Get number of particles per water molecule by inspecting the last residue in the topology
                resi_generator = modeller.topology.residues()
                resi_deque = deque(resi_generator, maxlen=1)
                last_resi = resi_deque.pop()
                nparticles_per_water = len([atom for atom in last_resi.atoms()])

                # Count final atoms.
                natoms_final = len(positions)
                nwaters = (natoms_final - natoms_initial) / nparticles_per_water
                if verbose: print("Solvated model contains %d waters" % nwaters)

                # Record waters.
                with open(nwaters_filename, 'w') as nwaters_file:
                    nwaters_file.write('%d\n' % nwaters)

            except Exception as e:
                reject_file_path = os.path.join(model_dir, 'solvation-rejected.txt')
                exception_text = '%r' % e
                trbk = traceback.format_exc()
                with open(reject_file_path, 'w') as reject_file:
                    reject_file.write(exception_text + '\n')
                    reject_file.write(trbk + '\n')

        if mpistate.rank == 0:
            project_metadata = ensembler.core.ProjectMetadata(project_stage='solvate_models', target_id=target.id)
            datestamp = ensembler.core.get_utcnow_formatted()
            target_timedelta = datetime.datetime.utcnow() - target_starttime

            metadata = {
                'target_id': target.id,
                'datestamp': datestamp,
                'template_seqid_cutoff': template_seqid_cutoff,
                'process_only_these_targets': process_only_these_targets,
                'process_only_these_templates': process_only_these_templates,
                'python_version': sys.version.split('|')[0].strip(),
                'python_full_version': ensembler.core.literal_str(sys.version),
                'ensembler_version': ensembler.version.short_version,
                'ensembler_commit': ensembler.version.git_revision,
                'biopython_version': Bio.__version__,
                'openmm_version': simtk.openmm.version.short_version,
                'openmm_commit': simtk.openmm.version.git_revision,
                'timing': ensembler.core.strf_timedelta(target_timedelta),
            }

            project_metadata.add_data(metadata)
            project_metadata.write()

        mpistate.comm.Barrier()

    mpistate.comm.Barrier()
    if mpistate.rank == 0:
        print('Done.')


def determine_nwaters(process_only_these_targets=None,
                      process_only_these_templates=None, template_seqid_cutoff=None,
                      verbose=False,
                      select_at_percentile=None):
    '''Determine distribution of nwaters, and select the value at a certain percentile.
    If not user-specified, the percentile is set to 100 if there are less than 10 templates, otherwise it is set to 68.
    '''

    # Run serially
    if mpistate.rank == 0:
        models_dir = os.path.abspath(ensembler.core.default_project_dirnames.models)

        targets, templates_resolved_seq = ensembler.core.get_targets_and_templates()

        if process_only_these_templates:
            selected_template_indices = [i for i, seq in enumerate(templates_resolved_seq) if seq.id in process_only_these_templates]
        else:
            selected_template_indices = range(len(templates_resolved_seq))

        for target in targets:

            # Process only specified targets if directed.
            if process_only_these_targets and (target.id not in process_only_these_targets): continue

            models_target_dir = os.path.join(models_dir, target.id)
            if not os.path.exists(models_target_dir): continue

            if template_seqid_cutoff:
                process_only_these_templates = ensembler.core.select_templates_by_seqid_cutoff(target.id, seqid_cutoff=template_seqid_cutoff)
                selected_template_indices = [i for i, seq in enumerate(templates_resolved_seq) if seq.id in process_only_these_templates]

            ntemplates_selected = len(selected_template_indices)

            if not select_at_percentile:
                select_at_percentile = 100 if ntemplates_selected < 10 else 68

            if verbose: print("Determining number of waters in each system from target '%s'..." % target.id)

            nwaters_list = []
            for template_index in range(ntemplates_selected):
                template = templates_resolved_seq[selected_template_indices[template_index]]
                if process_only_these_templates and template.id not in process_only_these_templates:
                    continue

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

            nwaters_array = np.array(nwaters_list)
            nwaters_array.sort()

            nwaters_list_filename = os.path.join(models_target_dir, 'nwaters-list.txt')
            with open(nwaters_list_filename, 'w') as nwaters_list_file:
                for nwaters in nwaters_array:
                    nwaters_list_file.write('%12d\n' % nwaters)

            # display statistics
            index_selected = int((len(nwaters_array) - 1) * (float(select_at_percentile) / 100.0))
            index68 = int((len(nwaters_array) - 1) * 0.68)
            index95 = int((len(nwaters_array) - 1) * 0.95)
            if len(nwaters_array) > 0:
                logger.info('Number of waters in solvated models (target: %s): min = %d, max = %d, '
                            'mean = %.1f, 68%% = %.0f, 95%% = %.0f, chosen_percentile (%d%%) = %.0f' %
                            (
                                target.id,
                                nwaters_array.min(),
                                nwaters_array.max(),
                                nwaters_array.mean(),
                                nwaters_array[index68],
                                nwaters_array[index95],
                                select_at_percentile,
                                nwaters_array[index_selected]
                            )
                            )

                filename = os.path.join(models_target_dir, 'nwaters-max.txt')
                with open(filename, 'w') as outfile:
                    outfile.write('%d\n' % nwaters_array.max())

                # Use 68th percentile.
                filename = os.path.join(models_target_dir, 'nwaters-use.txt')
                with open(filename, 'w') as outfile:
                    outfile.write('%d\n' % nwaters_array[index_selected])

            else:
                logger.info('No nwaters information found.')

            project_metadata = ensembler.core.ProjectMetadata(project_stage='determine_nwaters', target_id=target.id)

            datestamp = ensembler.core.get_utcnow_formatted()

            metadata = {
                'target_id': target.id,
                'datestamp': datestamp,
                'template_seqid_cutoff': template_seqid_cutoff,
                'select_at_percentile': select_at_percentile,
                'process_only_these_targets': process_only_these_targets,
                'process_only_these_templates': process_only_these_templates,
                'python_version': sys.version.split('|')[0].strip(),
                'python_full_version': ensembler.core.literal_str(sys.version),
                'ensembler_version': ensembler.version.short_version,
                'ensembler_commit': ensembler.version.git_revision,
                'biopython_version': Bio.__version__,
            }

            project_metadata.add_data(metadata)
            project_metadata.write()

        mpistate.comm.Barrier()

    mpistate.comm.Barrier()
    if mpistate.rank == 0:
        print('Done.')


def refine_explicit_md(
        openmm_platform=None, gpupn=1, process_only_these_targets=None,
        process_only_these_templates=None, template_seqid_cutoff=None,
        verbose=False, write_trajectory=False,
        include_disulfide_bonds=False,
        ff='amber99sbildn',
        water_model='tip3p',
        sim_length=100.0 * unit.picoseconds,
        timestep=2.0 * unit.femtoseconds, # timestep
        temperature=300.0 * unit.kelvin, # simulation temperature
        pressure=1.0 * unit.atmospheres, # simulation pressure
        collision_rate=20.0 / unit.picoseconds, # Langevin collision rate
        barostat_period=50,
        minimization_tolerance=10.0 * unit.kilojoules_per_mole / unit.nanometer,
        minimization_steps=20,
        nsteps_per_iteration=500,
        write_solvated_model=False,
        cpu_platform_threads=1,
        retry_failed_runs=False,
        serialize_at_start_of_each_sim=False):
    '''Run MD refinement in explicit solvent.

    MPI-enabled.
    '''
    gpuid = mpistate.rank % gpupn

    models_dir = os.path.abspath(ensembler.core.default_project_dirnames.models)

    targets, templates_resolved_seq = ensembler.core.get_targets_and_templates()

    if process_only_these_templates:
        selected_template_indices = [i for i, seq in enumerate(templates_resolved_seq) if seq.id in process_only_these_templates]
    else:
        selected_template_indices = range(len(templates_resolved_seq))

    if not openmm_platform:
        openmm_platform = auto_select_openmm_platform()

    if openmm_platform == 'CPU':
        platform_properties = {'CpuThreads': str(cpu_platform_threads)}
    else:
        platform_properties = {}

    ff_files = [ff+'.xml', water_model+'.xml']
    forcefield = app.ForceField(*ff_files)

    nonbondedMethod = app.PME

    kB = unit.MOLAR_GAS_CONSTANT_R
    kT = kB * temperature

    niterations = int((sim_length / timestep) / nsteps_per_iteration)

    def solvate_pdb(pdb, target_nwaters, water_model=water_model):
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

        There is no error checking to be sure that waters are not initially present in the system or that initially-present molecules are not deleted.

        """

        natoms_per_solvent = 3

        # Count initial atoms.
        natoms_initial = len(pdb.positions)
        if verbose: print("System initially has %d atoms (0 waters)" % (natoms_initial))

        # Solvate with zero padding to determine min number of waters and minimal unit cell dimensions.
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addSolvent(forcefield, model=water_model, padding=0.0*unit.angstroms)
        topology = modeller.getTopology()
        positions = modeller.getPositions()
        box_min = topology.getUnitCellDimensions()
        natoms_min = len(positions) # minimal number of atoms
        nwaters_min = (natoms_min - natoms_initial) / natoms_per_solvent # minimal number of waters
        volume_min = box_min[0] * box_min[1] * box_min[2]
        residues = [ r for r in topology.residues() ] # build a list of residues
        nresidues_min = len(residues) # number of residues
        if verbose: print("Minimally solvated system has %d atoms (%d waters)" % (natoms_min, nwaters_min))

        # Increase the box size by 10% and resolvate.
        scale = 1.1
        modeller = app.Modeller(pdb.topology, pdb.positions)
        topology = modeller.getTopology()
        topology.setUnitCellDimensions(box_min * scale)
        modeller.addSolvent(forcefield, model=water_model)
        positions = modeller.getPositions()
        box_enlarged = topology.getUnitCellDimensions()
        natoms_enlarged = len(positions) # minimal number of atoms
        nwaters_enlarged = (natoms_enlarged - natoms_initial) / natoms_per_solvent # minimal number of waters
        volume_enlarged = box_enlarged[0] * box_enlarged[1] * box_enlarged[2]
        density = (nwaters_enlarged - nwaters_min) / (volume_enlarged - volume_min)
        if verbose: print("Enlarged solvated system has %d atoms (%d waters) : density of %.3f waters / nm^3" % (natoms_enlarged, nwaters_enlarged, density / (1.0 / unit.nanometer**3)))

        # Aim for slightly more waters than target.
        over_target = False
        extra_nwaters = 100
        while not over_target:
            delta_volume = (target_nwaters + extra_nwaters - nwaters_min) / density
            scale = ((volume_min + delta_volume) / volume_min)**(1.0/3.0)
            if verbose: print("Final target of %d waters, so attempting box size %s to achieve %d waters..." % (target_nwaters, str(box_min * scale), target_nwaters + extra_nwaters))
            delta_volume = (target_nwaters + extra_nwaters - nwaters_min) / density
            modeller = app.Modeller(pdb.topology, pdb.positions)
            topology = modeller.getTopology()
            topology.setUnitCellDimensions(box_min * scale)
            modeller.addSolvent(forcefield, model=water_model)
            positions = modeller.getPositions()
            topology = modeller.getTopology()
            natoms = len(positions) # minimal number of atoms
            nwaters = (natoms - natoms_initial) / natoms_per_solvent # minimal number of waters
            if verbose: print("  actual %d waters" % nwaters)
            if (nwaters > target_nwaters):
                over_target = True
            else:
                extra_nwaters += 100

        # Delete waters to achieve target.
        ndelete = nwaters - target_nwaters
        if (ndelete > 0):
            if verbose: print("Will delete %d waters..." % ndelete)
            residues = [ r for r in topology.residues() ] # build a list of residues
            nresidues = len(residues)

            # Select a random subset to delete.
            indices = np.random.permutation(range(nresidues_min,nresidues))
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
        else:
            if write_solvated_model:
                # write solvated pdb file
                with open(os.path.join(model_dir, 'model-solvated.pdb'), 'w') as pdb_outfile:
                    app.PDBFile.writeHeader(topology, file=pdb_outfile)
                    app.PDBFile.writeFile(topology, positions, file=pdb_outfile)
                    app.PDBFile.writeFooter(topology, file=pdb_outfile)

        return [positions, topology]

    def simulate_explicit_md():
        # Set up Platform
        platform = openmm.Platform.getPlatformByName(openmm_platform)
        if 'CUDA_VISIBLE_DEVICES' not in os.environ:
            # Set GPU id.
            if openmm_platform == 'CUDA':
                platform.setPropertyDefaultValue('CudaDeviceIndex', '%d' % gpuid)
            elif openmm_platform == 'OpenCL':
                platform.setPropertyDefaultValue('OpenCLDeviceIndex', '%d' % gpuid)

        if verbose: print("Constructing System object...")
        system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, constraints=app.HBonds)
        if verbose: print("  system has %d atoms" % (system.getNumParticles()))

        # Add barostat.
        if verbose: print("Adding barostat...")
        barostat = openmm.MonteCarloBarostat(pressure, temperature, barostat_period)
        system.addForce(barostat)

        if verbose: print("Creating Context...")
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator, platform, platform_properties)
        context.setPositions(positions)

        if verbose: print("Minimizing structure...")
        openmm.LocalEnergyMinimizer.minimize(context, minimization_tolerance, minimization_steps)

        if write_trajectory:
            # Open trajectory for writing.
            if verbose: print("Opening trajectory for writing...")
            trajectory_filename = os.path.join(model_dir, 'explicit-trajectory.pdb.gz')
            trajectory_outfile = gzip.open(trajectory_filename, 'w')
            app.PDBFile.writeHeader(pdb.topology, file=trajectory_outfile)

        # Open energy trajectory for writing
        energy_filename = os.path.join(model_dir, 'explicit-energies.txt')
        energy_outfile = open(energy_filename, 'w')
        energy_outfile.write('# iteration | simulation time (ps) | potential_energy (kT) | kinetic_energy (kT) | volume (nm^3) | ns per day\n')

        if verbose: print("Running dynamics...")
        context.setVelocitiesToTemperature(temperature)
        import time
        initial_time = time.time()

        # TODO DEBUGGING
        if serialize_at_start_of_each_sim:
            with open(system_filename[: system_filename.index('.xml')]+'-start.xml', 'w') as system_file:
                system_file.write(openmm.XmlSerializer.serialize(system))
            with open(integrator_filename[: integrator_filename.index('.xml')]+'-start.xml', 'w') as integrator_file:
                integrator_file.write(openmm.XmlSerializer.serialize(integrator))
            state = context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
            with open(state_filename[: state_filename.index('.xml')]+'-start.xml', 'w') as state_file:
                state_file.write(openmm.XmlSerializer.serialize(state))

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
            if verbose: print("  %8.1f ps : potential %8.3f kT | kinetic %8.3f kT | volume %.3f nm^3 | %.3f ns/day | %.3f s remain" % (simulation_time / unit.picoseconds, potential_energy / kT, kinetic_energy / kT, volume_in_nm3, ns_per_day, remaining_time / unit.seconds))

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
        with gzip.open(pdb_filename, 'w') as pdb_outfile:
            app.PDBFile.writeHeader(topology, file=pdb_outfile)
            app.PDBFile.writeFile(topology, state.getPositions(), file=pdb_outfile)
            app.PDBFile.writeFooter(topology, file=pdb_outfile)

        # Serialize system
        if verbose: print("Serializing system...")
        with gzip.open(system_filename+'.gz', 'w') as system_file:
            system_file.write(openmm.XmlSerializer.serialize(system))

        # Serialize integrator.
        if verbose: print("Serializing integrator...")
        with gzip.open(integrator_filename+'.gz', 'w') as integrator_file:
            integrator_file.write(openmm.XmlSerializer.serialize(integrator))

        # Serialize state.
        if verbose: print("Serializing state...")
        state = context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
        with gzip.open(state_filename+'.gz', 'w') as state_file:
            state_file.write(openmm.XmlSerializer.serialize(state))


    for target in targets:
        if process_only_these_targets and (target.id not in process_only_these_targets):
            continue
        models_target_dir = os.path.join(models_dir, target.id)
        if mpistate.rank == 0:
            target_starttime = datetime.datetime.utcnow()
            if not os.path.exists(models_target_dir):
                continue

        mpistate.comm.Barrier()

        # Determine number of waters to use.
        nwaters_filename = os.path.join(models_target_dir, 'nwaters-use.txt')
        with open(nwaters_filename, 'r') as infile:
            line = infile.readline()
        nwaters = int(line)

        if template_seqid_cutoff:
            process_only_these_templates = ensembler.core.select_templates_by_seqid_cutoff(target.id, seqid_cutoff=template_seqid_cutoff)
            selected_template_indices = [i for i, seq in enumerate(templates_resolved_seq) if seq.id in process_only_these_templates]

        ntemplates_selected = len(selected_template_indices)

        for template_index in range(mpistate.rank, ntemplates_selected, mpistate.size):
            template = templates_resolved_seq[selected_template_indices[template_index]]

            model_dir = os.path.join(models_target_dir, template.id)
            if not os.path.exists(model_dir): continue

            # Only simulate models that are unique following filtering by clustering.
            unique_by_clustering = os.path.exists(os.path.join(model_dir, 'unique_by_clustering'))
            if not unique_by_clustering: continue

            # Pass if this simulation has already been run.
            log_filepath = os.path.join(model_dir, 'explicit-log.yaml')
            if os.path.exists(log_filepath):
                with open(log_filepath) as log_file:
                    try:
                        log_data = yaml.load(log_file, Loader=ensembler.core.YamlLoader)
                        if log_data.get('successful') is True:
                            continue
                        if log_data.get('finished') is True and (retry_failed_runs is False and log_data.get('successful') is False):
                            continue
                    except ScannerError as e:
                        trbk = traceback.format_exc()
                        warnings.warn(
                            '= WARNING start: template {0} MPI rank {1} hostname {2} gpuid {3} =\n{4}\n{5}\n= WARNING end: template {0} MPI rank {1} hostname {2} gpuid {3}'.format(
                                template.id, mpistate.rank, socket.gethostname(), gpuid, e, trbk
                            )
                        )

            # Check to make sure the initial model file is present.
            model_filename = os.path.join(model_dir, 'implicit-refined.pdb.gz')
            if not os.path.exists(model_filename):
                if verbose: print('model.pdb.gz not present: target %s template %s rank %d gpuid %d' % (target.id, template.id, mpistate.rank, gpuid))
                continue

            pdb_filename = os.path.join(model_dir, 'explicit-refined.pdb.gz')
            system_filename = os.path.join(model_dir, 'explicit-system.xml')
            integrator_filename = os.path.join(model_dir, 'explicit-integrator.xml')
            state_filename = os.path.join(model_dir, 'explicit-state.xml')

            print("-------------------------------------------------------------------------")
            print("Simulating %s => %s in explicit solvent for %.1f ps" % (target.id, template.id, niterations * nsteps_per_iteration * timestep / unit.picoseconds))
            print("-------------------------------------------------------------------------")

            # Open log file
            log_data = {
                'mpi_rank': mpistate.rank,
                'gpuid': gpuid if 'CUDA_VISIBLE_DEVICES' not in os.environ else os.environ['CUDA_VISIBLE_DEVICES'],
                'openmm_platform': openmm_platform,
                'sim_length': '%s' % sim_length,
                'finished': False,
                }
            log_file = ensembler.core.LogFile(log_filepath)
            log_file.log(new_log_data=log_data)

            try:
                start = datetime.datetime.utcnow()

                with gzip.open(model_filename) as model_file:
                    pdb = app.PDBFile(model_file)

                if not include_disulfide_bonds:
                    remove_disulfide_bonds_from_topology(pdb.topology)

                [positions, topology] = solvate_pdb(pdb, nwaters)

                simulate_explicit_md()

                timing = ensembler.core.strf_timedelta(datetime.datetime.utcnow() - start)
                log_data = {
                    'finished': True,
                    'timing': timing,
                    'successful': True,
                    }
                log_file.log(new_log_data=log_data)

            except Exception as e:
                trbk = traceback.format_exc()
                warnings.warn(
                    '= ERROR start: template {0} MPI rank {1} hostname {2} gpuid {3} =\n{4}\n{5}\n= ERROR end: template {0} MPI rank {1} hostname {2} gpuid {3}'.format(
                        template.id, mpistate.rank, socket.gethostname(), gpuid, e, trbk
                    )
                )
                timing = ensembler.core.strf_timedelta(datetime.datetime.utcnow() - start)
                log_data = {
                    'exception': e,
                    'traceback': ensembler.core.literal_str(trbk),
                    'timing': timing,
                    'finished': True,
                    'successful': False,
                    }
                log_file.log(new_log_data=log_data)

        if verbose:
            print('Finished template loop: rank %d' % mpistate.rank)

        mpistate.comm.Barrier()

        if mpistate.rank == 0:
            project_metadata = ensembler.core.ProjectMetadata(project_stage='refine_explicit_md', target_id=target.id)
            datestamp = ensembler.core.get_utcnow_formatted()
            nsuccessful_refinements = subprocess.check_output(['find', models_target_dir, '-name', 'explicit-refined.pdb.gz']).count('\n')
            target_timedelta = datetime.datetime.utcnow() - target_starttime

            metadata = {
                'target_id': target.id,
                'datestamp': datestamp,
                'template_seqid_cutoff': template_seqid_cutoff,
                'process_only_these_targets': process_only_these_targets,
                'process_only_these_templates': process_only_these_templates,
                'timing': ensembler.core.strf_timedelta(target_timedelta),
                'ff': ff,
                'water_model': water_model,
                'nsuccessful_refinements': nsuccessful_refinements,
                'python_version': sys.version.split('|')[0].strip(),
                'python_full_version': ensembler.core.literal_str(sys.version),
                'ensembler_version': ensembler.version.short_version,
                'ensembler_commit': ensembler.version.git_revision,
                'biopython_version': Bio.__version__,
                'openmm_version': simtk.openmm.version.short_version,
                'openmm_commit': simtk.openmm.version.git_revision
            }

            project_metadata.add_data(metadata)
            project_metadata.write()

        mpistate.comm.Barrier()

    mpistate.comm.Barrier()
    if mpistate.rank == 0:
        print('Done.')


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

