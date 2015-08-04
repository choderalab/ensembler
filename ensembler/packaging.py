import os
import subprocess
from ensembler.core import mpistate, logger, default_project_dirnames
from ensembler.core import get_targets_and_templates, select_templates_by_seqid_cutoff
from ensembler.utils import set_loglevel, read_file_contents_gz_or_not
from ensembler.refinement import auto_select_openmm_platform
import simtk.unit as unit
import simtk.openmm as openmm

fah_projects_dir = os.path.join(default_project_dirnames.packaged_models, 'fah-projects')


def package_for_fah(process_only_these_targets=None,
                    process_only_these_templates=None,
                    template_seqid_cutoff=None,
                    nclones=1, archive=False,
                    openmm_platform=None,
                    loglevel=None):
    """
    Create the input files and directory structure necessary to start a Folding@Home project.

    MPI-enabled.

    Parameters
    ----------
    archive : Bool
        A .tgz compressed archive will be created for each individual RUN directory.
    """
    set_loglevel(loglevel)

    if mpistate.rank == 0:
        if not os.path.exists(fah_projects_dir):
            os.mkdir(fah_projects_dir)
    mpistate.comm.Barrier()

    targets, templates_resolved_seq = get_targets_and_templates()

    if not openmm_platform:
        openmm_platform = auto_select_openmm_platform()

    for target in targets:
        if process_only_these_targets and (target.id not in process_only_these_targets):
            continue

        target_project_dir = os.path.join(fah_projects_dir, target.id)

        models_target_dir = os.path.join(default_project_dirnames.models, target.id)
        if not os.path.exists(models_target_dir):
            continue

        mpistate.comm.Barrier()

        sorted_valid_templates = []

        if mpistate.rank == 0:
            logger.info('-------------------------------------------------------------------------')
            logger.info('Building FAH OpenMM project for target {}'.format(target.id))
            logger.info('-------------------------------------------------------------------------')

            valid_templates = get_valid_templates_for_target(
                target,
                templates_resolved_seq,
                process_only_these_templates,
                template_seqid_cutoff
            )

            sorted_valid_templates = sort_valid_templates_by_seqid(
                target,
                valid_templates
            )

            create_target_project_dir(target)

        sorted_valid_templates = mpistate.comm.bcast(sorted_valid_templates, root=0)

        logger.debug("Building RUNs in parallel...")

        for run_index in range(mpistate.rank, len(sorted_valid_templates), mpistate.size):
            logger.info('-------------------------------------------------------------------------')
            logger.info('Building RUN for template {}'.format(sorted_valid_templates[run_index].id))
            logger.info('-------------------------------------------------------------------------')

            template = sorted_valid_templates[run_index]

            source_dir = os.path.join(models_target_dir, template.id)
            generate_fah_run(
                target_project_dir,
                template,
                source_dir,
                run_index,
                nclones,
                openmm_platform,
            )

            if archive:
                archive_fah_run(target, run_index)

    mpistate.comm.Barrier()
    if mpistate.rank == 0:
        print('Done.')


filenames_necessary_for_fah_packaging = [
    'unique_by_clustering',
    'sequence-identity.txt',
    'explicit-system.xml',
    'explicit-state.xml',
    'explicit-integrator.xml',
]


def get_valid_templates_for_target(target,
                                   templates_resolved_seq,
                                   process_only_these_templates,
                                   template_seqid_cutoff
                                   ):
    logger.debug("Building list of valid templates...")
    models_target_dir = os.path.join(default_project_dirnames.models, target.id)
    if template_seqid_cutoff:
        selected_templates = select_templates_by_seqid_cutoff(
            target.id, seqid_cutoff=template_seqid_cutoff
        )
    elif process_only_these_templates:
        selected_templates = [
            seq_obj for seq_obj in templates_resolved_seq
            if seq_obj.id in process_only_these_templates
        ]
    else:
        selected_templates = templates_resolved_seq

    valid_templates = []

    for template in selected_templates:
        # Check to make sure all files needed are present.
        for filename in filenames_necessary_for_fah_packaging:
            fullpath = os.path.join(models_target_dir, template.id, filename)
            if not (os.path.exists(fullpath) or os.path.exists(fullpath+'.gz')):
                continue
        valid_templates.append(template)

    logger.debug('{} valid unique initial starting conditions found'.format(len(valid_templates)))

    return valid_templates


def sort_valid_templates_by_seqid(target, valid_templates):
    logger.debug("Sorting templates in order of decreasing sequence identity...")
    models_target_dir = os.path.join(default_project_dirnames.models, target.id)

    seqids = []

    for template in valid_templates:
        seqids.append(get_seqid_for_model(models_target_dir, template))

    sorted_valid_templates_and_seqids = sorted(
        zip(valid_templates, seqids),
        reverse=True,
        key=lambda x: x[1]
    )

    sorted_valid_templates = zip(*sorted_valid_templates_and_seqids)[0]
    return sorted_valid_templates


def get_seqid_for_model(models_target_dir, template):
    seqid_filename = os.path.join(models_target_dir, template.id, 'sequence-identity.txt')
    with open(seqid_filename, 'r') as infile:
        seqid = float(infile.readline().strip())
    return seqid


def create_target_project_dir(target):
    target_project_dir = os.path.join(fah_projects_dir, target.id)
    if not os.path.exists(target_project_dir):
        os.makedirs(target_project_dir)


def generate_fah_run(target_project_dir,
                     template,
                     source_dir,
                     run_index,
                     nclones,
                     openmm_platform,
                     ):
    """
    Build Folding@Home RUN and CLONE subdirectories from (possibly compressed) OpenMM serialized XML files.

    ARGUMENTS

    run (int) - run index
    """
    logger.debug("Building RUN %d" % run_index)

    try:
        # Determine directory and pathnames.
        run_dir = os.path.join(target_project_dir, 'RUN%d' % run_index)
        run_template_id_filepath = os.path.join(run_dir, 'template.txt')
        run_seqid_filepath = os.path.join(run_dir, 'sequence-identity.txt')
        run_system_filepath = os.path.join(run_dir, 'system.xml')
        run_integrator_filepath = os.path.join(run_dir, 'integrator.xml')
        run_protein_structure_filepath = os.path.join(run_dir, 'protein.pdb')
        run_system_structure_filepath = os.path.join(run_dir, 'system.pdb')
        run_final_state_filepath = os.path.join(run_dir, 'state%d.xml' % (nclones - 1))
        source_seqid_filepath = os.path.join(source_dir, 'sequence-identity.txt')
        source_protein_structure_filepath = os.path.join(source_dir, 'implicit-refined.pdb.gz')
        source_system_structure_filepath = os.path.join(source_dir, 'explicit-refined.pdb.gz')
        source_openmm_system_filepath = os.path.join(source_dir, 'explicit-system.xml')
        source_openmm_state_filepath = os.path.join(source_dir, 'explicit-state.xml')

        # Return if this directory has already been set up.
        if os.path.exists(run_dir):
            if (
                    os.path.exists(run_template_id_filepath)
                    and os.path.exists(run_seqid_filepath)
                    and os.path.exists(run_system_filepath)
                    and os.path.exists(run_integrator_filepath)
                    and os.path.exists(run_protein_structure_filepath)
                    and os.path.exists(run_system_structure_filepath)
                    and os.path.exists(run_final_state_filepath)
                    ):
                return
        else:
            # Construct run directory if it does not exist.
            if not os.path.exists(run_dir):
                os.makedirs(run_dir)

        # Write template ID
        with open(run_template_id_filepath, 'w') as outfile:
            outfile.write(template.id + '\n')

        # Write the protein and system structure pdbs

        with open(run_protein_structure_filepath, 'w') as protein_structure_file:
            protein_structure_file.write(
                read_file_contents_gz_or_not(source_protein_structure_filepath)
            )

        with open(run_system_structure_filepath, 'w') as system_structure_file:
            system_structure_file.write(
                read_file_contents_gz_or_not(source_system_structure_filepath)
            )

        system = openmm.XmlSerializer.deserialize(
            read_file_contents_gz_or_not(source_openmm_system_filepath)
        )
        state = openmm.XmlSerializer.deserialize(
            read_file_contents_gz_or_not(source_openmm_state_filepath)
        )

        # Substitute default box vectors.
        box_vectors = state.getPeriodicBoxVectors()
        system.setDefaultPeriodicBoxVectors(*box_vectors)

        # Write sequence identity.
        with open(run_seqid_filepath, 'w') as run_seqid_file:
            run_seqid_file.write(read_file_contents_gz_or_not(source_seqid_filepath))

        # Integrator settings.
        constraint_tolerance = 1.0e-5
        timestep = 2.0 * unit.femtoseconds
        collision_rate = 1.0 / unit.picosecond
        temperature = 300.0 * unit.kelvin

        # Create new integrator to use.
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)

        # TODO: Make sure MonteCarloBarostat temperature matches set temperature.

        # Serialize System.
        with open(run_system_filepath, 'w') as run_system_file:
            run_system_file.write(openmm.XmlSerializer.serialize(system))

        # Serialize Integrator
        with open(run_integrator_filepath, 'w') as run_integrator_file:
            run_integrator_file.write(openmm.XmlSerializer.serialize(integrator))

        # Create Context so we can randomize velocities.
        platform = openmm.Platform.getPlatformByName(openmm_platform)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(state.getPositions())
        context.setVelocities(state.getVelocities())
        box_vectors = state.getPeriodicBoxVectors()
        context.setPeriodicBoxVectors(*box_vectors)

        # Create clones with different random initial velocities.
        for clone_index in range(nclones):
            state_filename = os.path.join(run_dir, 'state%d.xml' % clone_index)
            if os.path.exists(state_filename):
                continue
            context.setVelocitiesToTemperature(temperature)
            state = context.getState(
                getPositions=True,
                getVelocities=True,
                getForces=True,
                getEnergy=True,
                getParameters=True,
                enforcePeriodicBox=True
            )
            with open(state_filename, 'w') as state_file:
                state_file.write(openmm.XmlSerializer.serialize(state))

    except Exception as e:
        import traceback
        print(traceback.format_exc())
        print(str(e))

    return


def archive_fah_run(target, run_index):
    project_target_dir = os.path.join(fah_projects_dir, target.id)
    archive_filename = os.path.join(project_target_dir, 'RUN%d.tgz' % run_index)
    run_dir = os.path.join(project_target_dir, 'RUN%d' % run_index)
    subprocess.call(['tar', 'zcf', archive_filename, run_dir])


def package_for_transfer(process_only_these_targets=None):
    raise Exception('Not implemented yet.')
