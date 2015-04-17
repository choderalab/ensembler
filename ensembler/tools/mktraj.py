import os
import numpy as np
import pandas as pd
import mdtraj
import ensembler
from ensembler.core import logger, check_ensembler_modeling_stage_complete


def mktraj(targetid, ensembler_stage=None, traj_filepath=None, topol_filepath=None, models_data_filepath=None, process_only_these_templates=None):
    """Makes a trajectory for a given target, using mdtraj. The trajectory can be used with other
    software, e.g. for visualization with PyMOL or VMD.

    Parameters
    ----------
    targetid : str
        e.g. 'EGFR_HUMAN_D0'
    ensembler_stage : str
        The Ensembler stage from which to build models, e.g. 'build_models' results in a trajectory
        built from the 'model.pdb.gz' files output by the build_models command.
        options: build_models|refine_implicit_md|refine_explicit_md
        default: most advanced stage for which model files are available
    traj_filepath : str
        default: models/[targetid]/traj-[ensembler_stage].xtc
    topol_filepath : str
        default: models/[targetid]/traj-[ensembler_stage]-topol.pdb
    models_data_filepath :
        default: models/[targetid]/traj-[ensembler_stage]-data.csv
    process_only_these_templates : list of str

    Returns
    -------
    traj : mdtraj.Trajectory
    df : pandas.DataFrame
        models data (e.g. sequence identities):
    """
    try:
        ensembler.core.check_project_toplevel_dir()
        models_target_dir = os.path.join(ensembler.core.default_project_dirnames.models, targetid)

        logger.debug('Working on target %s' % targetid)

        if ensembler_stage is None:
            for stagename in ['refine_explicit_md', 'refine_implicit_md', 'build_models']:
                if check_ensembler_modeling_stage_complete(stagename, targetid):
                    ensembler_stage = stagename
                    break

        if ensembler_stage is None:
            raise Exception('Models have not yet been built for this Ensembler project.')

        if traj_filepath is None:
            traj_filepath = os.path.join(models_target_dir, 'traj-{0}.xtc'.format(ensembler_stage))
        if topol_filepath is None:
            topol_filepath = os.path.join(models_target_dir, 'traj-{0}-topol.pdb'.format(ensembler_stage))
        if models_data_filepath is None:
            models_data_filepath = os.path.join(models_target_dir, 'traj-{0}-data.csv'.format(ensembler_stage))

        if process_only_these_templates:
            templateids = process_only_these_templates
        else:
            dirs = os.walk(models_target_dir).next()[1]
            templateids = [dir for dir in dirs if '_D' in dir]

        model_filename = ensembler.core.model_filenames_by_ensembler_stage[ensembler_stage]
        valid_model_templateids = [templateid for templateid in templateids if os.path.exists(os.path.join(models_target_dir, templateid, model_filename))]
        valid_model_filepaths = [os.path.join(models_target_dir, templateid, model_filename) for templateid in valid_model_templateids]

        seqid_filepaths = [os.path.join(models_target_dir, templateid, 'sequence-identity.txt') for templateid in valid_model_templateids]
        seqids = [float(open(seqid_filepath).read().strip()) if os.path.exists(seqid_filepath) else None for seqid_filepath in seqid_filepaths]

        df = pd.DataFrame({
            'templateid': valid_model_templateids,
            'model_filepath': valid_model_filepaths,
            'seqid': seqids,
        })
        df.sort(columns='seqid', inplace=True, ascending=False)
        df.reset_index(drop=True, inplace=True)

        df.to_csv(models_data_filepath, columns=['templateid', 'seqid'])

        # construct traj
        traj = mdtraj.load_pdb(df.model_filepath[0])

        for model_filepath in df.model_filepath[1:]:
            traj += mdtraj.load_pdb(model_filepath)

        # superpose structured C-alphas
        dssp = mdtraj.compute_dssp(traj[0])[0]
        structured_resis_bool = (dssp == 'H') + (dssp == 'E')
        alpha_indices = traj.topology.select_atom_indices('alpha')
        structured_alpha_indices = np.array([alpha_indices[x] for x in range(traj.n_residues) if structured_resis_bool[x]])
        traj.superpose(reference=traj, frame=0, atom_indices=structured_alpha_indices)

        # write traj, and write first frame as pdb file
        traj[0].save(topol_filepath)
        traj.save(traj_filepath)
        return traj, df

    except:
        import ipdb; ipdb.set_trace()