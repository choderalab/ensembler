import os
import numpy as np
import pandas as pd
import mdtraj
import ensembler
from ensembler.core import logger, get_most_advanced_ensembler_modeling_stage


class MkTraj(object):
    def __init__(self, targetid, ensembler_stage=None, traj_filepath=None, topol_filepath=None,
           models_data_filepath=None, process_only_these_templates=None, loglevel=None,
           run_main=True):
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
        ensembler.utils.set_loglevel(loglevel)
        ensembler.core.check_project_toplevel_dir()
        self.models_target_dir = os.path.join(ensembler.core.default_project_dirnames.models, targetid)

        logger.debug('Working on target %s' % targetid)

        if ensembler_stage is None:
            self.ensembler_stage = get_most_advanced_ensembler_modeling_stage(targetid)
        else:
            self.ensembler_stage = ensembler_stage

        if traj_filepath is None:
            self.traj_filepath = os.path.join(
                self.models_target_dir, 'traj-{0}.xtc'.format(ensembler_stage)
            )
        else:
            self.traj_filepath = traj_filepath

        if topol_filepath is None:
            self.topol_filepath = os.path.join(
                self.models_target_dir, 'traj-{0}-topol.pdb'.format(ensembler_stage)
            )
        else:
            self.topol_filepath = topol_filepath

        if models_data_filepath is None:
            self.models_data_filepath = os.path.join(
                self.models_target_dir, 'traj-{0}-data.csv'.format(ensembler_stage)
            )
        else:
            self.models_data_filepath = models_data_filepath

        if process_only_these_templates:
            self.templateids = process_only_these_templates
        else:
            self.templateids = os.walk(self.models_target_dir).next()[1]

        if run_main:
            self._gen_df()
            self.df.to_csv(self.models_data_filepath, columns=['templateid', 'seqid'])
            self._construct_traj()
            self._superpose()
            self._write_traj()

    def _gen_df(self):
        model_filename = ensembler.core.model_filenames_by_ensembler_stage[self.ensembler_stage]
        valid_model_templateids = [
            templateid for templateid in self.templateids
            if os.path.exists(os.path.join(self.models_target_dir, templateid, model_filename))
        ]
        valid_model_filepaths = [
            os.path.join(self.models_target_dir, templateid, model_filename)
            for templateid in valid_model_templateids
        ]

        seqid_filepaths = [
            os.path.join(self.models_target_dir, templateid, 'sequence-identity.txt')
            for templateid in valid_model_templateids
        ]
        seqids = [
            float(open(seqid_filepath).read().strip()) if os.path.exists(seqid_filepath) else None
            for seqid_filepath in seqid_filepaths
        ]

        self.df = pd.DataFrame({
            'templateid': valid_model_templateids,
            'model_filepath': valid_model_filepaths,
            'seqid': seqids,
        })
        self.df.sort(columns='seqid', inplace=True, ascending=False)
        self.df.reset_index(drop=True, inplace=True)

    def _construct_traj(self):
        logger.debug('Working on model {0} ({1}/{2})'.format(self.df.templateid.iloc[0], 0, len(self.df.model_filepath)))
        traj = mdtraj.load_pdb(self.df.model_filepath[0])
        self.traj = traj

        for m, model_filepath in enumerate(self.df.model_filepath[1:]):
            logger.debug('Working on model {0} ({1}/{2})'.format(self.df.templateid.iloc[m+1], m+1, len(self.df.model_filepath)))
            traj = mdtraj.load_pdb(model_filepath)
            self.traj += traj

    def _superpose(self):
        """
        Superpose structured C-alphas
        """
        self.dssp = mdtraj.compute_dssp(self.traj[0])[0]
        structured_resis_bool = (self.dssp == 'H') + (self.dssp == 'E')
        alpha_indices = self.traj.topology.select_atom_indices('alpha')
        structured_alpha_indices = np.array([
            alpha_indices[x] for x in range(self.traj.n_residues) if structured_resis_bool[x]
        ])
        self.traj.superpose(reference=self.traj, frame=0, atom_indices=structured_alpha_indices)

    def _write_traj(self):
        """
        Write traj, and write first frame as pdb file
        """
        self.traj[0].save(self.topol_filepath)
        self.traj.save(self.traj_filepath)