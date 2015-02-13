import os
import yaml
import numpy as np
import mdtraj
import ensembler
from ensembler.core import logger


def mk_traj(targetid, traj_filepath=None, process_only_these_templates=None):
    models_target_dir = os.path.join(ensembler.core.default_project_dirnames.models, targetid)

    # target_meta_filepath = os.path.join(models_target_dir, 'meta.yaml')
    # with open(target_meta_filepath) as target_meta_file:
    #     target_meta = yaml.load(target_meta_file)
    #     if 'build_models' not in target_meta.keys():
    #         continue

    logger.debug('Working on target %s' % targetid)

    if traj_filepath is None:
        traj_filepath = os.path.join(models_target_dir, 'traj-models.h5')

    if process_only_these_templates:
        templateids = process_only_these_templates
    else:
        # get template names from sequence identities file
        seqids_filepath = os.path.join(models_target_dir, 'sequence-identities.txt')
        with open(seqids_filepath) as seqids_file:
            seqids = [x.split()[0] for x in seqids_file.readlines()]
        model_filepaths = [os.path.join(models_target_dir, seqid, 'model.pdb.gz') for seqid in seqids]
        templateids = []

    model_filepaths = [os.path.join(models_target_dir, templateid, 'model.pdb.gz') for templateid in templateids]
    valid_model_filepaths = [model_filepath for model_filepath in model_filepaths if os.path.exists(model_filepath)]

    # construct traj
    models = []
    traj = mdtraj.load_pdb(valid_model_filepaths[0])
    for model_filepath in valid_model_filepaths[1:]:
        traj = traj.join(mdtraj.load_pdb(model_filepath))

    # superpose structured C-alphas
    dssp = mdtraj.compute_dssp(traj[0])[0]
    structured_resis_bool = (dssp == 'H') + (dssp == 'E')
    alpha_indices = traj.topology.select_atom_indices('alpha')
    structured_alpha_indices = np.array([alpha_indices[x] for x in range(traj.n_residues) if structured_resis_bool[x]])
    traj.superpose(reference=traj, frame=0, atom_indices=structured_alpha_indices)

    traj.save(traj_filepath)