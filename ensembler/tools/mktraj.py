import os
import gzip
import numpy as np
import pandas as pd
import mdtraj
import ensembler
from ensembler.core import logger, get_most_advanced_ensembler_modeling_stage, default_project_dirnames, model_filenames_by_ensembler_stage, mpistate
from ensembler.refinement import remove_disulfide_bonds_from_topology, get_highest_seqid_existing_model


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
        self.models_target_dir = os.path.join(default_project_dirnames.models, targetid)

        logger.debug('Working on target %s' % targetid)

        if ensembler_stage is None:
            self.ensembler_stage = get_most_advanced_ensembler_modeling_stage(targetid)
        else:
            self.ensembler_stage = ensembler_stage

        if traj_filepath is None:
            self.traj_filepath = os.path.join(
                self.models_target_dir, 'traj-{0}.xtc'.format(self.ensembler_stage)
            )
        else:
            self.traj_filepath = traj_filepath

        if topol_filepath is None:
            self.topol_filepath = os.path.join(
                self.models_target_dir, 'traj-{0}-topol.pdb'.format(self.ensembler_stage)
            )
        else:
            self.topol_filepath = topol_filepath

        if models_data_filepath is None:
            self.models_data_filepath = os.path.join(
                self.models_target_dir, 'traj-{0}-data.csv'.format(self.ensembler_stage)
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

    def _gen_df(self, model_filename=None):
        if model_filename is None:
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
        remove_disulfide_bonds_from_topology(traj.topology)
        self.traj = traj

        for m, model_filepath in enumerate(self.df.model_filepath[1:]):
            logger.debug('Working on model {0} ({1}/{2})'.format(self.df.templateid.iloc[m+1], m+1, len(self.df.model_filepath)))
            traj = mdtraj.load_pdb(model_filepath)
            remove_disulfide_bonds_from_topology(traj.topology)
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


class MkTrajImplicitStart(MkTraj):
    def __init__(self, targetid, traj_filepath=None, topol_filepath=None,
           models_data_filepath=None, process_only_these_templates=None, loglevel=None,
           run_main=True):
        """Quick hack.
        """
        ensembler.utils.set_loglevel(loglevel)
        ensembler.core.check_project_toplevel_dir()
        self.models_target_dir = os.path.join(default_project_dirnames.models, targetid)

        logger.debug('Working on target %s' % targetid)

        self.ensembler_stage = 'implicit-start'
        self.model_filename = 'implicit-start.pdb.gz'

        if traj_filepath is None:
            self.traj_filepath = os.path.join(
                self.models_target_dir, 'traj-{0}.xtc'.format(self.ensembler_stage)
            )
        else:
            self.traj_filepath = traj_filepath

        if topol_filepath is None:
            self.topol_filepath = os.path.join(
                self.models_target_dir, 'traj-{0}-topol.pdb'.format(self.ensembler_stage)
            )
        else:
            self.topol_filepath = topol_filepath

        if models_data_filepath is None:
            self.models_data_filepath = os.path.join(
                self.models_target_dir, 'traj-{0}-data.csv'.format(self.ensembler_stage)
            )
        else:
            self.models_data_filepath = models_data_filepath

        if process_only_these_templates:
            self.templateids = process_only_these_templates
        else:
            self.templateids = os.walk(self.models_target_dir).next()[1]

        if run_main:
            self._gen_implicit_start_models()
            self._gen_df(model_filename=self.model_filename)
            self.df.to_csv(self.models_data_filepath, columns=['templateid', 'seqid'])
            self._construct_traj()
            self._superpose()
            self._write_traj()

    def _gen_implicit_start_models(
            self,
            ff='amber99sbildn.xml', implicit_water_model='amber99_obc.xml',
            ph=8.0):

        self.ph = ph
        from simtk.openmm import app

        valid_model_templateids = [
            templateid for templateid in self.templateids
            if os.path.exists(
                os.path.join(
                    self.models_target_dir, templateid,
                    ensembler.core.model_filenames_by_ensembler_stage['refine_implicit_md']
                )
            )
        ]

        gen_model_templateids = [
            templateid for templateid in valid_model_templateids
            if not os.path.exists(
                os.path.join(self.models_target_dir, templateid, self.model_filename)
            )
        ]

        # make reference model
        forcefield = app.ForceField(ff, implicit_water_model)
        reference_model_id = get_highest_seqid_existing_model(models_target_dir=self.models_target_dir)
        reference_model_path = os.path.join(self.models_target_dir, reference_model_id, model_filenames_by_ensembler_stage['build_models'])
        with gzip.open(reference_model_path) as reference_pdb_file:
            reference_pdb = app.PDBFile(reference_pdb_file)
        remove_disulfide_bonds_from_topology(reference_pdb.topology)
        reference_topology = reference_pdb.topology
        reference_modeller = app.Modeller(reference_pdb.topology, reference_pdb.positions)
        reference_variants = reference_modeller.addHydrogens(forcefield, pH=self.ph)

        for template_index in range(mpistate.rank, len(gen_model_templateids), mpistate.size):
            templateid = gen_model_templateids[template_index]
            logger.debug('Generating implicit-start model for {0}'.format(templateid))

            try:
                input_model_filepath = os.path.join(self.models_target_dir, templateid, model_filenames_by_ensembler_stage['build_models'])
                output_model_filepath = os.path.join(self.models_target_dir, templateid, self.model_filename)

                with gzip.open(input_model_filepath) as pdb_file:
                    pdb = app.PDBFile(pdb_file)

                remove_disulfide_bonds_from_topology(pdb.topology)
                modeller = app.Modeller(reference_topology, pdb.positions)
                modeller.addHydrogens(forcefield, pH=self.ph, variants=reference_variants)
                topology = modeller.getTopology()
                positions = modeller.getPositions()

                with gzip.open(output_model_filepath, 'w') as output_model_file:
                    app.PDBFile.writeHeader(topology, file=output_model_file)
                    app.PDBFile.writeFile(topology, positions, file=output_model_file)
                    app.PDBFile.writeFooter(topology, file=output_model_file)

            except Exception as e:
                continue
                # import ipdb; ipdb.set_trace()