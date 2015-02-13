import os
import yaml
import numpy as np
import pandas as pd
import mdtraj
import ensembler
from ensembler.core import logger


def mk_traj(targetid, traj_filepath=None, topol_filepath=None, models_data_filepath=None, process_only_these_templates=None):
    """Makes a trajectory for a given target, using mdtraj. The trajectory can be used with other
    software, e.g. for visualization with PyMOL or VMD.

    Parameters
    ----------
    targetid : str
        e.g. 'EGFR_HUMAN_D0'
    traj_filepath : str
        default: models/[targetid]/modelstraj.xtc
    topol_filepath : str
        default: models/[targetid]/modelstraj-topol.pdb
    models_data_filepath :
        default: models/[targetid]/modelstraj-data.csv
    process_only_these_templates : list of str

    Returns
    -------
    traj : mdtraj.Trajectory
    df : pandas.DataFrame
        models data (e.g. sequence identities):
    """
    models_target_dir = os.path.join(ensembler.core.default_project_dirnames.models, targetid)

    logger.debug('Working on target %s' % targetid)

    if traj_filepath is None:
        traj_filepath = os.path.join(models_target_dir, 'modelstraj.xtc')
    if topol_filepath is None:
        topol_filepath = os.path.join(models_target_dir, 'modelstraj-topol.pdb')
    if models_data_filepath is None:
        models_data_filepath = os.path.join(models_target_dir, 'modelstraj-data.csv')

    if process_only_these_templates:
        templateids = process_only_these_templates
    else:
        dirs = os.walk(models_target_dir).next()[1]
        templateids = [dir for dir in dirs if '_D' in dir]

    valid_model_templateids = [templateid for templateid in templateids if os.path.exists(os.path.join(models_target_dir, templateid, 'model.pdb.gz'))]
    valid_model_filepaths = [os.path.join(models_target_dir, templateid, 'model.pdb.gz') for templateid in valid_model_templateids]

    # get additional data
    seqid_filepaths = [os.path.join(models_target_dir, templateid, 'sequence-identity.txt') for templateid in valid_model_templateids]
    seqids = [open(seqid_filepath).read().strip() if os.path.exists(seqid_filepath) else None for seqid_filepath in seqid_filepaths]

    df = pd.DataFrame({
        'templateid': valid_model_templateids,
        'seqid': seqids,
    })
    df.to_csv(models_data_filepath)

    # construct traj
    traj = mdtraj.load_pdb(valid_model_filepaths[0])
    for model_filepath in valid_model_filepaths[1:]:
        traj = traj.join(mdtraj.load_pdb(model_filepath))

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