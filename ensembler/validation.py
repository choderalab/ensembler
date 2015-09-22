import os
import tempfile
import shutil
import yaml
import mdtraj
from subprocess import Popen, PIPE
from ensembler.core import get_most_advanced_ensembler_modeling_stage, default_project_dirnames
from ensembler.core import model_filenames_by_ensembler_stage, get_valid_model_ids, mpistate
from ensembler.core import YamlDumper, YamlLoader, logger
from ensembler.utils import notify_when_done, set_loglevel

# includes types
molprobity_oneline_analysis_colnames = [
    ('#pdbFileName', None),
    ('x-H_type', None),
    ('chains', int),
    ('residues', int),
    ('nucacids', int),
    ('resolution', float),
    ('rvalue', float),
    ('rfree', float),
    ('clashscore', float),
    ('clashscoreB<40', float),
    ('minresol', float),
    ('maxresol', float),
    ('n_samples', int),
    ('pct_rank', int),
    ('pct_rank40', int),
    ('cbeta>0.25', int),
    ('numCbeta', int),
    ('rota<1%', int),
    ('numRota', int),
    ('ramaOutlier', int),
    ('ramaAllowed', int),
    ('ramaFavored', int),
    ('numRama', int),
    ('numbadbonds', int),
    ('numbonds', int),
    ('pct_badbonds', float),
    ('pct_resbadbonds', float),
    ('numbadangles', int),
    ('numangles', int),
    ('pct_badangles', float),
    ('pct_resbadangles', float),
    ('MolProbityScore', float),
    ('Mol_pct_rank', int),
]

molprobity_oneline_analysis_colnames_to_output = [
    'MolProbityScore',
    'clashscore',
    'numRota',
    'rota<1%',
    'numRama',
    'ramaOutlier',
    'ramaFavored',
    'cbeta>0.25',
    'pct_badbonds',
    'pct_badangles',
]


@notify_when_done
def molprobity_validation_multiple_targets(targetids, modeling_stage=None, loglevel=None):
    """
Calculate model quality using MolProbity ``oneline-analysis`` command.

For each target, this function outputs a text file named
``models/[targetid]/validation_scores_sorted-[method]-[ensembler_stage]`` which contains a list of
targetids sorted by validation score. This can be used by the subsequent ``package_models`` command
to filter out models below a specified quality threshold.

Typically, this should be run after models have been refined to the desired extent (e.g. after
implicit or explicit MD refinement)

More detailed validation results are written to the individual model directories.

MPI-enabled.

    Parameters
    ----------
    targetids: list of str or str
    modeling_stage: str
        {None|build_models|refine_implicit_md|refine_explicit_md}
        Default: None (automatically selects most advanced stage)
    """
    set_loglevel(loglevel)
    if type(targetids) is str:
        targetids = [targetids]
    for targetid in targetids:
        logger.info('Working on target {}'.format(targetid))
        molprobity_validation(targetid=targetid, ensembler_stage=modeling_stage, loglevel=loglevel)


def molprobity_validation(targetid, ensembler_stage=None, loglevel=None):
    set_loglevel(loglevel)
    valid_model_ids = []
    if mpistate.rank == 0:
        if ensembler_stage is None:
            ensembler_stage = get_most_advanced_ensembler_modeling_stage(targetid)
        valid_model_ids = get_valid_model_ids(ensembler_stage, targetid)
    if ensembler_stage is None:
        ensembler_stage = mpistate.comm.bcast(ensembler_stage, root=0)
    valid_model_ids = mpistate.comm.bcast(valid_model_ids, root=0)
    nvalid_model_ids = len(valid_model_ids)
    model_structure_filename = model_filenames_by_ensembler_stage[ensembler_stage]

    models_target_dir = os.path.join(default_project_dirnames.models, targetid)
    molprobity_results_filepath = os.path.join(
        models_target_dir, 'validation_scores_sorted-molprobity-{}'.format(ensembler_stage)
    )

    molprobity_scores_sublist = []
    for model_index in range(mpistate.rank, nvalid_model_ids, mpistate.size):
        model_id = valid_model_ids[model_index]

        logger.debug('MPI process {} working on model {}'.format(mpistate.rank, model_id))

        molprobity_score = run_molprobity_oneline_analysis_and_write_results(
            targetid,
            model_id,
            ensembler_stage,
            model_structure_filename=model_structure_filename,
            models_target_dir=models_target_dir,
        )

        molprobity_scores_sublist.append((model_id, molprobity_score))

    molprobity_scores_gathered_list = mpistate.comm.gather(molprobity_scores_sublist, root=0)
    molprobity_scores_list_of_tuples = [item for sublist in molprobity_scores_gathered_list for item in sublist]
    molprobity_scores_list_of_tuples = mpistate.comm.bcast(molprobity_scores_list_of_tuples, root=0)
    molprobity_scores_sorted = sorted(molprobity_scores_list_of_tuples, key=lambda x: x[1])
    write_molprobity_scores_list(molprobity_scores_sorted, molprobity_results_filepath)


def run_molprobity_oneline_analysis_and_write_results(targetid,
                                                      model_id,
                                                      ensembler_stage,
                                                      model_structure_filename=None,
                                                      models_target_dir=None,
                                                      check_for_existing_results=True,
                                                      ):
    if model_structure_filename is None:
        model_structure_filename = model_filenames_by_ensembler_stage[ensembler_stage]
    if models_target_dir is None:
        models_target_dir = os.path.join(default_project_dirnames.models, targetid)

    results_output_filepath = os.path.join(
        models_target_dir, model_id, 'molprobity-{}.yaml'.format(ensembler_stage)
    )

    if check_for_existing_results:
        if os.path.exists(results_output_filepath):
            with open(results_output_filepath) as results_output_file:
                prev_results = yaml.load(stream=results_output_file, Loader=YamlLoader)
            prev_molprobity_score = prev_results.get('MolProbityScore')
            if prev_molprobity_score is not None:
                logger.debug(
                    'Existing MolProbity score of {} found for model {}'.format(
                        prev_molprobity_score, model_id
                    )
                )
                return prev_molprobity_score

    molprobity_results = run_molprobity_oneline_analysis(
        targetid, model_id, model_structure_filename
    )
    if molprobity_results is None:
        logger.debug('MolProbity returned no results for model {}'.format(model_id))
        return None

    logger.debug('MolProbity score of {} calculated for model {}'.format(molprobity_results.get('MolProbityScore'), model_id))

    molprobity_score = molprobity_results.get('MolProbityScore')
    if molprobity_score is not None:
        write_molprobity_results_for_target(
            molprobity_results, models_target_dir, model_id, ensembler_stage
        )
    return molprobity_score


def run_molprobity_oneline_analysis(targetid, model_id, model_structure_filename, tmp_model_dir=None):
    """
    Runs oneline_analysis for a single model in a temp dir, and cleans up after.
    """
    if tmp_model_dir is None:
        tmp_model_dir = tempfile.mkdtemp()

    try:
        source_path = os.path.join(
            default_project_dirnames.models,
            targetid,
            model_id,
            model_structure_filename
        )
        dest_path = os.path.join(
            tmp_model_dir,
            model_id + '.pdb'
        )

        source_model_traj = mdtraj.load_pdb(source_path)
        protein_only_traj = source_model_traj.atom_slice(
            source_model_traj.top.select('protein')
        )
        protein_only_traj.save_pdb(dest_path)

        stdout, stderr = molprobity_oneline_analysis_cmd(tmp_model_dir)
        output_text = '\n'.join([stdout, stderr])
        molprobity_results = parse_molprobity_oneline_analysis_output(output_text)
        molprobity_model_results = molprobity_results.get(model_id)

    finally:
        shutil.rmtree(tmp_model_dir)

    return molprobity_model_results


def molprobity_oneline_analysis_cmd(dir_path):
    p = Popen(
        [
            'oneline-analysis',
            dir_path
        ],
        stdout=PIPE,
        stderr=PIPE,
    )
    stdout, stderr = p.communicate()
    return stdout, stderr


def parse_molprobity_oneline_analysis_output(output_text):
    results_lines = []
    for line in output_text.splitlines():
        if len(line) == 0 or line[0] == '#':
            continue
        ncolons = line.count(':')
        if ncolons == 32:
            results_lines.append(line)

    molprobity_results = {}
    for results_line in results_lines:
        results_line_split = results_line.split(':')
        filename = results_line_split[0]
        targetid = filename[: filename.find('.pdb')]

        target_results = {}
        for c, coltuple in enumerate(molprobity_oneline_analysis_colnames):
            colname, coltype = coltuple
            value = results_line_split[c]
            try:
                if coltype is not None:
                    value = coltype(value)
            except (ValueError, TypeError):
                pass
            target_results[colname] = value

        molprobity_results[targetid] = target_results
    return molprobity_results


def write_molprobity_results_for_target(molprobity_model_results, models_target_dir, model_id, ensembler_stage):
    output_dict = {
        colname: molprobity_model_results[colname] for colname in molprobity_oneline_analysis_colnames_to_output
    }

    results_output_filepath = os.path.join(
        models_target_dir, model_id, 'molprobity-{}.yaml'.format(ensembler_stage)
    )
    with open(results_output_filepath, 'w') as results_output_file:
        yaml.dump(output_dict, stream=results_output_file, default_flow_style=False, Dumper=YamlDumper)


def write_molprobity_scores_list(molprobity_scores_sorted, molprobity_results_filepath):
    output_text = '\n'.join(
        ['{} {}'.format(score_tuple[0], score_tuple[1]) for score_tuple in molprobity_scores_sorted]
    )
    with open(molprobity_results_filepath, 'w') as molprobity_results_file:
        molprobity_results_file.write(output_text)
