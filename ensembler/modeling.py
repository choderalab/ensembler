import sys
import os
import re
import datetime
import shutil
import gzip
import glob
from collections import namedtuple
import tempfile
import traceback
import Bio.SeqUtils
import simtk.openmm
import yaml
import warnings
import ensembler
import ensembler.version
import Bio
import Bio.SeqIO
import Bio.pairwise2
import Bio.SubsMat.MatrixInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import mdtraj
import msmbuilder.cluster
from ensembler.core import get_targets_and_templates, get_templates_full_seq
from ensembler.core import mpistate, logger
try:
    import modeller
    import modeller.automodel
except ImportError:
    pass
try:
    import subprocess32 as subprocess
    loopmodel_subprocess_kwargs = {'timeout': 10800}   # 3 hour timeout - used for loopmodel call
except ImportError:
    warnings.warn('subprocess32 module not available. Falling back to subprocess module, without timeout functionality.')
    import subprocess
    loopmodel_subprocess_kwargs = {}


TargetSetupData = namedtuple(
    'TargetSetupData',
    ['target_starttime', 'models_target_dir']
)


class LoopmodelOutput:
    def __init__(self, output_text=None, loopmodel_exception=None, exception=None, trbk=None, successful=False, no_missing_residues=False):
        self.output_text = output_text
        self.exception = exception
        self.loopmodel_exception = loopmodel_exception
        self.traceback = trbk
        self.successful = successful
        self.no_missing_residues = no_missing_residues


@ensembler.utils.notify_when_done
def model_template_loops(process_only_these_templates=None, overwrite_structures=False, loglevel=None):
    """
    Use Rosetta loopmodel to model missing loops in template structures.
    Completed templates are stored in templates/structures-modeled-loops

    :param process_only_these_templates: list of str
    :param loglevel: str
    :return:
    """
    ensembler.utils.set_loglevel(loglevel)
    targets, templates_resolved_seq = ensembler.core.get_targets_and_templates()
    templates_full_seq = get_templates_full_seq()
    missing_residues_list = pdbfix_templates(templates_full_seq, process_only_these_templates=process_only_these_templates, overwrite_structures=overwrite_structures)
    loopmodel_templates(templates_resolved_seq, missing_residues_list, process_only_these_templates=process_only_these_templates, overwrite_structures=overwrite_structures)


def pdbfix_templates(templates_full_seq, process_only_these_templates=None, overwrite_structures=False):
    """
    Parameters
    ----------
    templates_full_seq: list of BioPython SeqRecord
        full UniProt sequence for span of the template (including unresolved residues)
    process_only_these_templates: list of str
    overwrite_structures: bool
    Returns
    -------
    missing_residues_list: list of list of OpenMM Residue
    """
    missing_residues_sublist = []
    ntemplates = len(templates_full_seq)
    for template_index in range(mpistate.rank, ntemplates, mpistate.size):
        template_full_seq = templates_full_seq[template_index]
        if process_only_these_templates and template_full_seq.id not in process_only_these_templates:
            missing_residues_sublist.append(None)
            continue
        missing_residues_sublist.append(pdbfix_template(template_full_seq, overwrite_structures=overwrite_structures))

    missing_residues_gathered = mpistate.comm.gather(missing_residues_sublist, root=0)

    missing_residues_list = []
    if mpistate.rank == 0:
        missing_residues_list = [None] * ntemplates
        for template_index in range(ntemplates):
            missing_residues_list[template_index] = missing_residues_gathered[template_index % mpistate.size][template_index // mpistate.size]

    missing_residues_list = mpistate.comm.bcast(missing_residues_list, root=0)

    return missing_residues_list


def pdbfix_template(template_full_seq, overwrite_structures=False):
    """
    Parameters
    ----------
    template_full_seq: BioPython SeqRecord
        full UniProt sequence for span of the template (including unresolved residues)
    overwrite_structures: bool
    Returns
    -------
    fixer.missingResidues
    """
    try:
        template_pdbfixed_filepath = os.path.join(
            ensembler.core.default_project_dirnames.templates_structures_modeled_loops,
            template_full_seq.id + '-pdbfixed.pdb'
        )
        seq_pdbfixed_filepath = os.path.join(
            ensembler.core.default_project_dirnames.templates_structures_modeled_loops,
            template_full_seq.id + '-pdbfixed.fasta'
        )
        import pdbfixer
        import simtk.openmm.app
        template_filepath = os.path.join(
            ensembler.core.default_project_dirnames.templates_structures_resolved,
            template_full_seq.id + '.pdb'
        )
        fixer = pdbfixer.PDBFixer(filename=template_filepath)
        chainid = fixer.structureChains[0].chain_id
        seq_obj = simtk.openmm.app.internal.pdbstructure.Sequence(chainid)
        for r in template_full_seq.seq:
            resi3 = Bio.SeqUtils.seq3(r).upper()
            seq_obj.residues.append(resi3)
        fixer.structure.sequences.append(seq_obj)
        fixer.findMissingResidues()
        remove_missing_residues_at_termini(fixer, len_full_seq=len(template_full_seq.seq))
        if not overwrite_structures and os.path.exists(template_pdbfixed_filepath):
            return fixer.missingResidues
        fixer.findMissingAtoms()
        (newTopology, newPositions, newAtoms, existingAtomMap) = fixer._addAtomsToTopology(True, True)
        fixer.topology = newTopology
        fixer.positions = newPositions
        with open(template_pdbfixed_filepath, 'w') as template_pdbfixed_file:
            simtk.openmm.app.PDBFile.writeFile(
                fixer.topology, fixer.positions, file=template_pdbfixed_file
            )

        # Write sequence to file
        seq_pdbfixed = ''.join([Bio.SeqUtils.seq1(r.name) for r in fixer.topology.residues()])
        seq_record_pdbfixed = SeqRecord(Seq(seq_pdbfixed), id=template_full_seq.id, description=template_full_seq.id)
        Bio.SeqIO.write([seq_record_pdbfixed], seq_pdbfixed_filepath, 'fasta')

        return fixer.missingResidues
    except (KeyboardInterrupt, ImportError):
        raise
    except Exception as e:
        trbk = traceback.format_exc()
        log_filepath = os.path.abspath(os.path.join(
            ensembler.core.default_project_dirnames.templates_structures_modeled_loops,
            template_full_seq.id + '-pdbfixer-log.yaml'
        ))
        logfile = ensembler.core.LogFile(log_filepath)
        logfile.log({
            'templateid': str(template_full_seq.id),
            'exception': e,
            'traceback': ensembler.core.literal_str(trbk),
            'mpi_rank': mpistate.rank,
        })
        logger.error(
            'MPI rank %d pdbfixer error for template %s - see logfile' %
            (mpistate.rank, template_full_seq.id)
        )


def remove_missing_residues_at_termini(fixer, len_full_seq):
    # remove C-terminal missing residues
    if len(fixer.missingResidues) == 0:
        return None
    sorted_missing_residues_keys = sorted(fixer.missingResidues, key=lambda x: x[1])
    last_missing_residues_key = sorted_missing_residues_keys[-1]
    last_missing_residues_start_index = last_missing_residues_key[1]
    last_missing_residues = fixer.missingResidues[last_missing_residues_key]
    nmissing_residues_up_to_last = sum([len(fixer.missingResidues[key]) for key in sorted_missing_residues_keys[:-1]])

    if last_missing_residues_start_index + nmissing_residues_up_to_last + len(last_missing_residues) == len_full_seq:
        fixer.missingResidues.pop(last_missing_residues_key)

    # remove N-terminal missing residues
    fixer.missingResidues.pop((0, 0), None)


def loopmodel_templates(templates, missing_residues, process_only_these_templates=None, overwrite_structures=False):
    """
    Parameters
    ----------
    templates:  list of BioPython SeqRecord
        only the id is used
    missing_residues: list of list of OpenMM Residue
    process_only_these_templates: bool
    overwrite_structures: bool
    """
    for template_index in range(mpistate.rank, len(templates), mpistate.size):
        template = templates[template_index]
        if process_only_these_templates and template.id not in process_only_these_templates:
            continue
        if mpistate.size > 1:
            logger.info('MPI rank %d modeling missing loops for template %s' % (mpistate.rank, template.id))
        else:
            logger.info('Modeling missing loops for template %s' % template.id)
        loopmodel_template(template, missing_residues[template_index], overwrite_structures=overwrite_structures)


def loopmodel_template(template, missing_residues, overwrite_structures=False):
    template_filepath = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.id + '-pdbfixed.pdb'))
    output_pdb_filepath = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.id + '.pdb'))
    loop_filepath = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.id + '.loop'))
    output_score_filepath = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.id + '-loopmodel-score.sc'))
    log_filepath = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.id + '-loopmodel-log.yaml'))
    if not overwrite_structures:
        if os.path.exists(log_filepath):
            return
    logfile = ensembler.core.LogFile(log_filepath)
    write_loop_file(template, missing_residues)
    starttime = datetime.datetime.utcnow()
    if len(missing_residues) == 0:
        loopmodel_output = LoopmodelOutput(successful=True, no_missing_residues=True)
    else:
        loopmodel_output = run_loopmodel(template_filepath, loop_filepath, output_pdb_filepath, output_score_filepath)
    if not loopmodel_output.successful:
        logger.error('MPI rank %d Loopmodel error for template %s - see logfile' % (mpistate.rank, template.id))
    timedelta = datetime.datetime.utcnow() - starttime
    logfile.log({
        'templateid': str(template.id),
        'no_missing_residues': loopmodel_output.no_missing_residues,
        'loopmodel_output': loopmodel_output.output_text,
        'mpi_rank': mpistate.rank,
        'successful': loopmodel_output.successful,
        'exception': loopmodel_output.exception,
        'loopmodel_exception': loopmodel_output.loopmodel_exception,
        'traceback': loopmodel_output.traceback,
        'timing': ensembler.core.strf_timedelta(timedelta),
        })


def write_loop_file(template, missing_residues):
    loop_file_text = ''
    loop_residues_added = 0
    loop_residues_data = [(key[1], len(residues)) for key, residues in missing_residues.iteritems()]
    loop_residues_data = sorted(loop_residues_data, key=lambda x: x[0])
    for loop_residue_data in loop_residues_data:
        residue_number, nresidues = loop_residue_data
        loop_begin = residue_number + loop_residues_added   # 1-based, one residue before the loop
        loop_end = residue_number + nresidues + loop_residues_added + 1   # 1-based, one residue after the loop
        loop_residues_added += nresidues
        # Note that missing residues at termini (which cannot be modeled by Rosetta loopmodel) have already been removed from the PDBFixer.missingResidues dictionary
        loop_file_text += 'LOOP {0} {1} - - 1\n'.format(loop_begin, loop_end)
    loop_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.id + '.loop')
    with open(loop_filepath, 'w') as loop_file:
        loop_file.write(loop_file_text)


def run_loopmodel(input_template_pdb_filepath, loop_filepath, output_pdb_filepath, output_score_filepath, loopmodel_executable_filepath=None, nmodels_to_build=1):
    if loopmodel_executable_filepath is None:
        loopmodel_executable_filepath = ensembler.core.find_loopmodel_executable()

    temp_dir = tempfile.mkdtemp()
    temp_template_filepath = os.path.join(temp_dir, 'template.pdb')
    temp_loop_filepath = os.path.join(temp_dir, 'template.loop')
    temp_output_model_filepath = os.path.join(temp_dir, 'template_0001.pdb')
    temp_output_score_filepath = os.path.join(temp_dir, 'score.sc')
    minirosetta_database_path = os.environ.get('MINIROSETTA_DATABASE')
    shutil.copy(input_template_pdb_filepath, temp_template_filepath)
    shutil.copy(loop_filepath, temp_loop_filepath)
    try:
        output_text = subprocess.check_output(
            [
                loopmodel_executable_filepath,
                '-database', minirosetta_database_path,
                '-in::file::s', temp_template_filepath,
                '-loops:loop_file', temp_loop_filepath,
                '-out:path:all', temp_dir,
                '-loops:remodel', 'perturb_kic',
                '-loops:refine', 'refine_kic',
                '-ex1',
                '-ex2',
                '-nstruct', '%d' % nmodels_to_build,
                '-loops:max_kic_build_attempts', '100',
                '-in:file:fullatom',
                '-overwrite',
                ],
            stderr=subprocess.STDOUT,
            **loopmodel_subprocess_kwargs
            )
        if os.path.exists(temp_output_model_filepath):
            shutil.copy(temp_output_model_filepath, output_pdb_filepath)
            shutil.copy(temp_output_score_filepath, output_score_filepath)
            shutil.rmtree(temp_dir)
            return LoopmodelOutput(output_text=output_text, successful=True)
        else:
            shutil.rmtree(temp_dir)
            return LoopmodelOutput(output_text=output_text, successful=False)
    except KeyboardInterrupt:
        shutil.rmtree(temp_dir)
        raise
    except subprocess.CalledProcessError as e:
        shutil.rmtree(temp_dir)
        return LoopmodelOutput(loopmodel_exception=e.output, trbk=traceback.format_exc(), successful=False)
    except subprocess.TimeoutExpired as e:
        shutil.rmtree(temp_dir)
        return LoopmodelOutput(output_text=e.output, exception=e, trbk=traceback.format_exc(), successful=False)
    except Exception as e:
        shutil.rmtree(temp_dir)
        return LoopmodelOutput(output_text=output_text, exception=e, trbk=traceback.format_exc(), successful=False)


def check_loopmodel_complete_and_successful(template):
    output_pdb_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.id + '.pdb')
    log_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.id + '-loopmodel-log.yaml')
    if os.path.exists(log_filepath) and os.path.exists(output_pdb_filepath):
        with open(log_filepath) as log_file:
            log_data = yaml.load(log_file, Loader=ensembler.core.YamlLoader)
            if log_data.get('successful') == True:
                return True
    else:
        return False


@ensembler.utils.notify_when_done
def align_targets_and_templates(process_only_these_targets=None, process_only_these_templates=None, loglevel=None):
    """
    Conducts pairwise alignments of target sequences against template sequences.
    Stores Modeller-compatible 'alignment.pir' files in each model directory,
    and also outputs a table of model IDs, sorted by sequence identity.

    :param process_only_these_targets:
    :param process_only_these_templates:
    :param loglevel:
    :return:
    """
    ensembler.utils.set_loglevel(loglevel)
    targets, templates_resolved_seq = ensembler.core.get_targets_and_templates()
    ntemplates = len(templates_resolved_seq)
    nselected_templates = len(process_only_these_templates) if process_only_these_templates else ntemplates
    for target in targets:
        if process_only_these_targets and target.id not in process_only_these_targets: continue

        if mpistate.rank == 0:
            logger.info('Working on target %s...' % target.id)

        models_target_dir = os.path.join(ensembler.core.default_project_dirnames.models, target.id)
        ensembler.utils.create_dir(models_target_dir)

        seq_identity_data_sublist = []

        for template_index in range(mpistate.rank, ntemplates, mpistate.size):
            template_id = templates_resolved_seq[template_index].id
            if os.path.exists(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template_id + '.pdb')):
                remodeled_seq_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template_id + '-pdbfixed.fasta')
                template = list(Bio.SeqIO.parse(remodeled_seq_filepath, 'fasta'))[0]
            else:
                template = templates_resolved_seq[template_index]

            if process_only_these_templates and template_id not in process_only_these_templates: continue

            model_dir = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.models, target.id, template_id))
            ensembler.utils.create_dir(model_dir)
            aln = align_target_template(target, template)
            aln_filepath = os.path.join(model_dir, 'alignment.pir')
            write_modeller_pir_aln_file(aln, target, template, pir_aln_filepath=aln_filepath)
            seq_identity_data_sublist.append({
                'templateid': template_id,
                'seq_identity': calculate_seq_identity(aln),
            })

        seq_identity_data_gathered = mpistate.comm.gather(seq_identity_data_sublist, root=0)

        seq_identity_data = []
        if mpistate.rank == 0:
            seq_identity_data = [None] * nselected_templates
            for i in range(nselected_templates):
                seq_identity_data[i] = seq_identity_data_gathered[i % mpistate.size][i // mpistate.size]

        seq_identity_data = mpistate.comm.bcast(seq_identity_data, root=0)

        seq_identity_data = sorted(seq_identity_data, key=lambda x: x['seq_identity'], reverse=True)
        write_sorted_seq_identities(target, seq_identity_data)


def align_target_template(target, template, gap_open=-10, gap_extend=-0.5):
    """
    :param target: BioPython SeqRecord
    :param template: BioPython SeqRecord
    :param gap_open: float or int
    :param gap_extend: float or int
    :return: alignment
    """
    matrix = Bio.SubsMat.MatrixInfo.gonnet
    aln = Bio.pairwise2.align.globalds(str(target.seq), str(template.seq), matrix, gap_open, gap_extend)
    return aln


def calculate_seq_identity(aln):
    len_shorter_seq = min([len(aln[0][0].replace('-', '')), len(aln[0][1].replace('-', ''))])
    seq_id = 0
    for r in range(len(aln[0][0])):
        if aln[0][0][r] == aln[0][1][r]:
            seq_id += 1
    seq_id = 100 * float(seq_id) / float(len_shorter_seq)
    return seq_id


def write_sorted_seq_identities(target, seq_identity_data):
    seq_identity_file_str = ''
    for seq_identity_dict in seq_identity_data:
        seq_identity_file_str += '%-30s %.1f\n' % (seq_identity_dict['templateid'], seq_identity_dict['seq_identity'])
    seq_identity_filepath = os.path.join(ensembler.core.default_project_dirnames.models, target.id, 'sequence-identities.txt')
    with open(seq_identity_filepath, 'w') as seq_identity_file:
        seq_identity_file.write(seq_identity_file_str)


@ensembler.utils.notify_when_done
def build_models(process_only_these_targets=None, process_only_these_templates=None,
                 template_seqid_cutoff=None, write_modeller_restraints_file=False, loglevel=None):
    """Uses the build_model method to build homology models for a given set of
    targets and templates.

    MPI-enabled.
    """
    # Note that this code uses an os.chdir call to switch into a temp directory before running Modeller.
    # This is because Modeller writes various output files in the current directory, and there is NO WAY
    # to define where these files are written, other than to chdir beforehand. If running this routine
    # in parallel, it is likely that occasional exceptions will occur, due to concurrent processes
    # making os.chdir calls.
    ensembler.utils.set_loglevel(loglevel)
    targets, templates_resolved_seq = get_targets_and_templates()

    if process_only_these_templates:
        selected_template_indices = [i for i, seq in enumerate(templates_resolved_seq) if seq.id in process_only_these_templates]
    else:
        selected_template_indices = range(len(templates_resolved_seq))

    for target in targets:
        if process_only_these_targets and target.id not in process_only_these_targets: continue
        target_setup_data = build_models_target_setup(target)

        if template_seqid_cutoff:
            process_only_these_templates = ensembler.core.select_templates_by_seqid_cutoff(target.id, seqid_cutoff=template_seqid_cutoff)
            selected_template_indices = [i for i, seq in enumerate(templates_resolved_seq) if seq.id in process_only_these_templates]

        ntemplates_selected = len(selected_template_indices)

        for template_index in range(mpistate.rank, ntemplates_selected, mpistate.size):
            template_resolved_seq = templates_resolved_seq[selected_template_indices[template_index]]
            if process_only_these_templates and template_resolved_seq.id not in process_only_these_templates: continue
            build_model(target, template_resolved_seq, target_setup_data,
                        write_modeller_restraints_file=write_modeller_restraints_file,
                        loglevel=loglevel)
        write_build_models_metadata(target, target_setup_data, process_only_these_targets,
                                    process_only_these_templates, template_seqid_cutoff,
                                    write_modeller_restraints_file)


def build_model(target, template_resolved_seq, target_setup_data,
                write_modeller_restraints_file=False, loglevel=None):
    """Uses Modeller to build a homology model for a given target and
    template.

    Will not run Modeller if the output files already exist.

    Parameters
    ----------
    target : BioPython SeqRecord
    template_resolved_seq : BioPython SeqRecord
        Must be a corresponding .pdb template file with the same ID in the
        templates/structures directory.
    template_resolved_seq : BioPython SeqRecord
        Must be a corresponding .pdb template file with the same ID in the
        templates/structures directory.
    target_setup_data : TargetSetupData obj
    write_modeller_restraints_file : bool
        Write file containing restraints used by Modeller - note that this file can be relatively
        large, e.g. ~300KB per model for a protein kinase domain target.
    loglevel : bool
    """
    ensembler.utils.set_loglevel(loglevel)

    template_structure_dir = os.path.abspath(
        ensembler.core.default_project_dirnames.templates_structures_modeled_loops
    )

    if os.path.exists(os.path.join(template_structure_dir, template_resolved_seq.id + '.pdb')):
        remodeled_seq_filepath = os.path.join(
            ensembler.core.default_project_dirnames.templates_structures_modeled_loops,
            template_resolved_seq.id + '-pdbfixed.fasta'
        )
        template = list(Bio.SeqIO.parse(remodeled_seq_filepath, 'fasta'))[0]
    else:
        template = template_resolved_seq
        template_structure_dir = os.path.abspath(
            ensembler.core.default_project_dirnames.templates_structures_resolved
        )

    model_dir = os.path.abspath(os.path.join(target_setup_data.models_target_dir, template.id))
    if not os.path.exists(model_dir):
        ensembler.utils.create_dir(model_dir)
    model_pdbfilepath = os.path.abspath(os.path.join(model_dir, 'model.pdb.gz'))
    modeling_log_filepath = os.path.abspath(os.path.join(model_dir, 'modeling-log.yaml'))

    check_model_pdbfilepath_ends_in_pdbgz(model_pdbfilepath)
    model_pdbfilepath_uncompressed = model_pdbfilepath[:-3]

    if check_all_model_files_present(model_dir):
        logger.debug(
            "Output files already exist for target '%s' // template '%s'; files were not overwritten." %
            (target.id, template.id)
        )
        return

    logger.info(
        '-------------------------------------------------------------------------\n'
        'Modelling "%s" => "%s"\n'
        '-------------------------------------------------------------------------'
        % (target.id, template.id)
    )

    # aln = align_target_template(target, template)
    aln_filepath = os.path.abspath(os.path.join(model_dir, 'alignment.pir'))
    # write_modeller_pir_aln_file(aln, target, template, pir_aln_filepath=aln_filepath)
    log_file = init_build_model_logfile(modeling_log_filepath)

    with ensembler.utils.enter_temp_dir():
        try:
            start = datetime.datetime.utcnow()
            shutil.copy(aln_filepath, 'alignment.pir')
            run_modeller(target, template, model_dir, model_pdbfilepath,
                         model_pdbfilepath_uncompressed, template_structure_dir,
                         write_modeller_restraints_file=write_modeller_restraints_file)
            if os.path.getsize(model_pdbfilepath) < 1:
                raise Exception('Output PDB file is empty.')

            end_successful_build_model_logfile(log_file, start)

        except Exception as e:
            end_exception_build_model_logfile(e, log_file)


def get_modeller_version():
    """Hacky attempt to get Modeller version by regex searching the installation directory or README file.
    """
    modeller_version = get_modeller_version_from_install_path(modeller)
    if modeller_version is not None:
        return modeller_version

    modeller_version = get_modeller_version_from_readme(modeller)
    if modeller_version is not None:
        return modeller_version


def get_modeller_version_from_install_path(modeller_module):
    regex = re.compile('/modeller-[0-9.]{2,6}/')
    match = re.search(regex, modeller_module.__file__)
    if match is not None:
        version = match.group()[10:-1]
        return version


def get_modeller_version_from_readme(modeller_module):
    readme_file_path = os.path.join(os.path.dirname(modeller_module.__file__), '..', '..', 'README')
    if os.path.exists(readme_file_path):
        with open(readme_file_path) as readme_file:
            # try first 10 lines
            # example desired line:
            #      MODELLER 9.11, 2012/08/29, r8834
            for i in range(10):
                line = readme_file.readline().strip()
                regex = re.compile('MODELLER [0-9.]{2,6}')
                match = re.search(regex, line)
                if match is not None:
                    version = match.group()[9:]
                    return version


def build_models_target_setup(target):
    target_setup_data = None
    if mpistate.rank == 0:
        models_target_dir = os.path.join(ensembler.core.default_project_dirnames.models, target.id)
        target_starttime = datetime.datetime.utcnow()
        logger.info(
            '=========================================================================\n'
            'Working on target "%s"\n'
            '========================================================================='
            % target.id
        )
        target_setup_data = TargetSetupData(
            target_starttime=target_starttime,
            models_target_dir=models_target_dir
        )
    target_setup_data = mpistate.comm.bcast(target_setup_data, root=0)
    return target_setup_data


def gen_build_models_metadata(target, target_setup_data, process_only_these_targets,
                              process_only_these_templates, template_seqid_cutoff,
                              write_modeller_restraints_file):
    """
    Generate build_models metadata for a given target.
    :param target: BioPython SeqRecord
    :param target_setup_data:
    :return: metadata: dict
    """
    datestamp = ensembler.core.get_utcnow_formatted()
    nsuccessful_models = subprocess.check_output(['find', target_setup_data.models_target_dir, '-name', 'model.pdb.gz']).count('\n')
    target_timedelta = datetime.datetime.utcnow() - target_setup_data.target_starttime
    modeller_version = get_modeller_version()
    metadata = {
        'target_id': target.id,
        'write_modeller_restraints_file': write_modeller_restraints_file,
        'template_seqid_cutoff': template_seqid_cutoff,
        'datestamp': datestamp,
        'timing': ensembler.core.strf_timedelta(target_timedelta),
        'nsuccessful_models': nsuccessful_models,
        'process_only_these_targets': process_only_these_targets,
        'process_only_these_templates': process_only_these_templates,
        'python_version': sys.version.split('|')[0].strip(),
        'python_full_version': ensembler.core.literal_str(sys.version),
        'ensembler_version': ensembler.version.short_version,
        'ensembler_commit': ensembler.version.git_revision,
        'modeller_version': modeller_version if modeller_version is not None else '',
        'biopython_version': Bio.__version__
    }
    return metadata


def check_model_pdbfilepath_ends_in_pdbgz(model_pdbfilepath):
    if model_pdbfilepath[-7:] != '.pdb.gz':
        raise Exception('model_pdbfilepath (%s) must end in .pdb.gz' % model_pdbfilepath)


def check_all_model_files_present(model_dir):
    seqid_filepath = os.path.abspath(os.path.join(model_dir, 'sequence-identity.txt'))
    model_pdbfilepath = os.path.abspath(os.path.join(model_dir, 'model.pdb.gz'))
    aln_filepath = os.path.abspath(os.path.join(model_dir, 'alignment.pir'))
    files_to_check = [model_pdbfilepath, seqid_filepath, aln_filepath]
    files_present = [os.path.exists(filename) for filename in files_to_check]
    return all(files_present)


def init_build_model_logfile(modeling_log_filepath):
    log_data = {
        'mpi_rank': mpistate.rank,
        'complete': False,
    }
    log_filepath = modeling_log_filepath
    log_file = ensembler.core.LogFile(log_filepath)
    log_file.log(new_log_data=log_data)
    return log_file


def write_modeller_pir_aln_file(aln, target, template, pir_aln_filepath='alignment.pir'):
    contents = "Target-template alignment\n"
    contents += ">P1;%s\n" % target.id
    contents += "sequence:%s:FIRST:@:LAST :@:::-1.00:-1.00\n" % target.id
    contents += aln[0][0] + '*\n'
    contents += ">P1;%s\n" % template.id
    contents += "structureX:%s:FIRST:@:LAST : :undefined:undefined:-1.00:-1.00\n" % template.id
    contents += aln[0][1] + '*\n'
    with open(pir_aln_filepath, 'w') as outfile:
        outfile.write(contents)


def run_modeller(target, template, model_dir, model_pdbfilepath, model_pdbfilepath_uncompressed,
                 template_structure_dir, aln_filepath='alignment.pir',
                 write_modeller_restraints_file=False):
    modeller.log.none()
    env = modeller.environ()
    env.io.atom_files_directory = [template_structure_dir]
    a = modeller.automodel.allhmodel(
        env,
        alnfile=aln_filepath,
        knowns=template.id,
        sequence=target.id
    )
    a.make()  # do homology modeling

    save_modeller_output_files(target, model_dir, a, env, model_pdbfilepath,
                               model_pdbfilepath_uncompressed,
                               write_modeller_restraints_file=write_modeller_restraints_file)


def save_modeller_output_files(target, model_dir, a, env, model_pdbfilepath,
                               model_pdbfilepath_uncompressed,
                               write_modeller_restraints_file=False):
    # save PDB file
    # Note that the uncompressed pdb file needs to be kept until after the clustering step has completed
    tmp_model_pdbfilepath = a.outputs[0]['name']
    target_model = modeller.model(env, file=tmp_model_pdbfilepath)
    target_model.write(file=model_pdbfilepath_uncompressed)
    with open(model_pdbfilepath_uncompressed) as model_pdbfile:
        with gzip.open(model_pdbfilepath, 'w') as model_pdbfilegz:
            model_pdbfilegz.write(model_pdbfile.read())

    # Write sequence identity.
    seqid_filepath = os.path.abspath(os.path.join(model_dir, 'sequence-identity.txt'))
    with open(seqid_filepath, 'w') as seqid_file:
        seqid_file.write('%.1f\n' % target_model.seq_id)

    # Copy restraints.
    if write_modeller_restraints_file:
        restraint_filepath = os.path.abspath(os.path.join(model_dir, 'restraints.rsr.gz'))
        with open('%s.rsr' % target.id, 'r') as rsrfile:
            with gzip.open(restraint_filepath, 'wb') as rsrgzfile:
                rsrgzfile.write(rsrfile.read())


def end_successful_build_model_logfile(log_file, start):
    end = datetime.datetime.utcnow()
    timing = ensembler.core.strf_timedelta(end - start)
    log_data = {
        'complete': True,
        'timing': timing,
    }
    log_file.log(new_log_data=log_data)


def end_exception_build_model_logfile(e, log_file):
    trbk = traceback.format_exc()
    log_data = {
        'exception': e,
        'traceback': ensembler.core.literal_str(trbk),
    }
    log_file.log(new_log_data=log_data)


@ensembler.utils.mpirank0only_and_end_with_barrier
def write_build_models_metadata(target, target_setup_data, process_only_these_targets,
                                process_only_these_templates, template_seqid_cutoff,
                                write_modeller_restraints_file):
    project_metadata = ensembler.core.ProjectMetadata(project_stage='build_models', target_id=target.id)
    metadata = gen_build_models_metadata(target, target_setup_data, process_only_these_targets,
                                         process_only_these_templates, template_seqid_cutoff,
                                         write_modeller_restraints_file)
    project_metadata.add_data(metadata)
    project_metadata.write()


@ensembler.utils.mpirank0only_and_end_with_barrier
@ensembler.utils.notify_when_done
def cluster_models(process_only_these_targets=None, cutoff=0.06, loglevel=None):
    """Cluster models based on RMSD, and filter out non-unique models as
    determined by a given cutoff.

    Parameters
    ----------

    cutoff : float
        Minimum distance cutoff for RMSD clustering (nm)

    Runs serially.
    """
    # TODO refactor
    ensembler.utils.set_loglevel(loglevel)
    targets, templates_resolved_seq = get_targets_and_templates()
    templates = templates_resolved_seq

    for target in targets:
        if process_only_these_targets and (target.id not in process_only_these_targets): continue

        models_target_dir = os.path.join(ensembler.core.default_project_dirnames.models, target.id)
        if not os.path.exists(models_target_dir): continue

        # =============================
        # Construct a mdtraj trajectory containing all models
        # =============================

        starttime = datetime.datetime.utcnow()

        logger.debug('Building a list of valid models...')

        model_pdbfilenames_compressed = {
            template.id: os.path.join(models_target_dir, template.id, 'model.pdb.gz') for template in templates
        }
        model_pdbfilenames_uncompressed = {
            template.id: os.path.join(models_target_dir, template.id, 'model.pdb') for template in templates
        }
        valid_templateids = [
            templateid for templateid in model_pdbfilenames_compressed
            if os.path.exists(model_pdbfilenames_compressed[templateid])
        ]

        # Write uncompressed model.pdb files from model.pdb.gz if necessary
        for templateid in valid_templateids:
            if not os.path.exists(model_pdbfilenames_uncompressed[templateid]) or os.path.getsize(model_pdbfilenames_uncompressed[templateid]) == 0:
                with gzip.open(model_pdbfilenames_compressed[templateid]) as model_pdbfile_compressed:
                    with open(model_pdbfilenames_uncompressed[templateid], 'w') as model_pdbfile:
                        model_pdbfile.write(model_pdbfile_compressed.read())

        logger.info('Constructing a trajectory containing all valid models...')

        if len(valid_templateids) == 0:
            logger.info('No models found for target {0}.'.format(target.id))
            continue

        valid_model_pdbfilenames_uncompressed = [
            model_pdbfilenames_uncompressed[templateid] for templateid in valid_templateids
        ]

        traj = mdtraj.load(valid_model_pdbfilenames_uncompressed)

        # =============================
        # Clustering
        # =============================

        logger.info('Conducting RMSD-based clustering...')

        # Remove any existing unique_by_clustering files
        for f in glob.glob(models_target_dir+'/*_PK_*/unique_by_clustering'):
            os.unlink(f)

        CAatoms = [a.index for a in traj.topology.atoms if a.name == 'CA']
        unique_templateids = models_regular_spatial_clustering(
            valid_templateids, traj, atom_indices=CAatoms, cutoff=cutoff
        )
        write_unique_by_clustering_files(unique_templateids, models_target_dir)

        with open(os.path.join(models_target_dir, 'unique-models.txt'), 'w') as uniques_file:
            for u in unique_templateids:
                uniques_file.write(u+'\n')
            logger.info(
                '%d unique models (from original set of %d) using cutoff of %.3f nm' %
                        (len(unique_templateids), len(valid_templateids), cutoff)
            )

        for template in templates:
            model_dir = os.path.join(models_target_dir, template.id)
            model_pdbfilename = os.path.join(model_dir, 'model.pdb')
            if os.path.exists(model_pdbfilename):
                os.remove(model_pdbfilename)

        # ========
        # Metadata
        # ========

        project_metadata = ensembler.core.ProjectMetadata(
            project_stage='cluster_models', target_id=target.id
        )
        datestamp = ensembler.core.get_utcnow_formatted()

        timedelta = datetime.datetime.utcnow() - starttime

        metadata = {
            'target_id': target.id,
            'datestamp': datestamp,
            'nunique_models': len(unique_templateids),
            'python_version': sys.version.split('|')[0].strip(),
            'python_full_version': ensembler.core.literal_str(sys.version),
            'ensembler_version': ensembler.version.short_version,
            'ensembler_commit': ensembler.version.git_revision,
            'biopython_version': Bio.__version__,
            'mdtraj_version': mdtraj.version.short_version,
            'mdtraj_commit': mdtraj.version.git_revision,
            'timing': ensembler.core.strf_timedelta(timedelta),
        }

        project_metadata.add_data(metadata)
        project_metadata.write()


def models_regular_spatial_clustering(templateids, traj, atom_indices=None, cutoff=0.06):
    """
    Use MSMBuilder to perform RMSD-based regular spatial clustering on a set of models.

    Parameters
    ----------
    templateids: list of str
    traj: mdtraj.Trajectory
    atom_indices: np.array
    cutoff: float
        Minimum distance cutoff for RMSD clustering (nm)
    """
    if atom_indices:
        reduced_traj = traj.atom_slice(atom_indices)
    else:
        reduced_traj = traj

    cluster = msmbuilder.cluster.RegularSpatial(cutoff, metric='rmsd')
    cluster_labels = cluster.fit_predict([reduced_traj])[0]
    unique_templateids = list(set([templateids[t] for t in cluster_labels]))
    return unique_templateids


def write_unique_by_clustering_files(unique_templateids, models_target_dir):
    for templateid in unique_templateids:
        unique_filename = os.path.join(models_target_dir, templateid, 'unique_by_clustering')
        with open(unique_filename, 'w') as unique_file:
            pass


def _deprecated_models_regular_spatial_clustering(templateids, traj, atom_indices=None, cutoff=0.06):
    """
    Superseded by models_regular_spatial_clustering
    """
    mdtraj_rmsd_args = {}
    if atom_indices:
        mdtraj_rmsd_args['atom_indices'] = atom_indices

    unique_templateids = []
    min_rmsd = []
    # Iterate through models
    for (t, templateid) in enumerate(templateids):
        # Add the first templateid to the list of uniques
        if t==0:
            unique_templateids.append(templateid)
            continue

        # Calculate rmsds of models up to t against the model t.
        rmsds = mdtraj.rmsd(traj[0:t], traj[t], parallel=False, **mdtraj_rmsd_args)
        min_rmsd.append(min(rmsds))

        # If any rmsd is less than cutoff, discard; otherwise add to list of uniques
        if min_rmsd[-1] < cutoff:
            continue
        else:
            unique_templateids.append(templateid)

    return unique_templateids