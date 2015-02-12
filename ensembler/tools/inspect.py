from datetime import datetime
import ensembler
import os
import pandas as pd
import yaml
import mdtraj
import gzip
from ensembler.core import logger
import warnings


class LoopmodelLogs(object):
    def __init__(self, project_dir='.'):
        self.project_dir = project_dir
        self.df = self.parse_loopmodel_logs()

    def parse_loopmodel_logs(self):
        templateid = []
        mpi_rank = []
        hostname = []
        successful = []
        no_missing_residues = []
        exception = []
        loopmodel_exception = []
        loopmodel_output = []
        traceback = []
        datestamp = []
        timing_str = []
        timing_total_seconds = []

        structures_modeled_loops_dir = os.path.join(self.project_dir, os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops))
        for filename in os.listdir(structures_modeled_loops_dir):
            if len(filename) >= 5:
                if filename[-5:] == '.yaml':
                    filepath = os.path.join(structures_modeled_loops_dir, filename)
                    with open(filepath) as logfile:
                        pass
                        log = yaml.load(logfile, Loader=ensembler.core.YamlLoader)

                        templateid.append(log.get('templateid'))
                        mpi_rank.append(int(log.get('mpi_rank')))
                        hostname.append(log.get('hostname'))
                        successful.append(log.get('successful'))
                        no_missing_residues.append(log.get('no_missing_residues'))
                        exception.append(log.get('exception'))
                        loopmodel_exception.append(log.get('loopmodel_exception'))
                        loopmodel_output.append(log.get('loopmodel_output'))
                        traceback.append(log.get('traceback'))
                        ds = datetime.strptime(log.get('datestamp'), ensembler.core.datestamp_format_string)
                        datestamp.append(ds)
                        timing_str.append(log.get('timing'))
                        timing_hours, timing_minutes, timing_seconds = map(int, timing_str[-1].split(':'))
                        timing_total_seconds.append(timing_hours*3600 + timing_minutes*60 + timing_seconds)

        df = pd.DataFrame(
            {
                'templateid': templateid,
                'mpi_rank': mpi_rank,
                'hostname': hostname,
                'successful': successful,
                'no_missing_residues': no_missing_residues,
                'exception': exception,
                'loopmodel_exception': loopmodel_exception,
                'loopmodel_output': loopmodel_output,
                'traceback': traceback,
                'datestamp': datestamp,
                'timing_str': timing_str,
                'timing_total_seconds': timing_total_seconds,
            }
        )
        return df

    def add_nmissing_resis_data(self):
        nmissing_resis_df = []
        for templateid in self.df.templateid:
            loopfile_path = os.path.join(self.project_dir, 'templates/structures-modeled-loops', templateid + '.loop')
            with open(loopfile_path) as loopfile:
                lines = loopfile.readlines()
                nmissing_resis = 0
                for line in lines:
                    loop_start = int(line[4:8])
                    loop_end = int(line[8:12])
                    nmissing_resis += loop_end - loop_start + 1
            nmissing_resis_df.append({
                'templateid': templateid,
                'nmissing_resis': nmissing_resis,
            })

        nmissing_resis_df = pd.DataFrame(nmissing_resis_df)

        self.df = self.df.merge(nmissing_resis_df, on='templateid')

    def to_hdf(self, ofilepath):
        self.df.to_hdf(ofilepath, 'df')

    def to_pickle(self, ofilepath):
        self.df.to_pickle(ofilepath)


class ModelSimilarities(object):
    def __init__(self, targetid, project_dir='.', log_level=None):
        ensembler.utils.loglevel_setter(logger, log_level)
        self.targetid = targetid
        self.project_dir = project_dir
        self._get_templateids_and_model_filepaths()
        self._get_unique_models()
        self._mk_traj()

    def _get_templateids_and_model_filepaths(self):
        model_dir = os.path.join(ensembler.core.default_project_dirnames.models, self.targetid)

        root, dirnames, filenames = os.walk(model_dir).next()

        templateids = [dirname for dirname in dirnames if '_D' in dirname]
        template_dirpaths = []
        has_model = []
        model_filepaths = []
        for templateid in templateids:
            template_dirpaths.append(os.path.join(root, templateid))
            model_filepath = os.path.join(root, templateid, 'model.pdb.gz')
            if os.path.exists(model_filepath):
                model_filepaths.append(model_filepath)
                has_model.append(True)
            else:
                has_model.append(False)

        self.templateids = templateids
        self.template_dirpaths = template_dirpaths
        self.model_filepaths = model_filepaths
        self.df = pd.DataFrame(
            {'templateid': templateids,
             'has_model': has_model}
        )

    def _get_unique_models(self):
        unique_models = []
        for template_dirpath in self.template_dirpaths:
            unique_path = os.path.join(template_dirpath, 'unique_by_clustering')
            if os.path.exists(unique_path):
                unique_models.append(True)
            else:
                unique_models.append(False)
        self.df['unique_by_clustering'] = unique_models

    def _mk_traj(self):
        with ensembler.utils.mk_temp_dir() as tmpdir:
            model_filepaths = []
            for m, model_filepath_gz in enumerate(self.model_filepaths):
                logger.debug('Unzipping model {0}/{1}'.format(m, len(self.model_filepaths)))
                with gzip.open(model_filepath_gz) as model_file:
                    model_filepath = os.path.join(tmpdir, '{0}.pdb'.format(m))
                    model_filepaths.append(model_filepath)
                    model_text = model_file.read()
                with open(model_filepath, 'w') as model_file:
                    model_file.write(model_text)

            self.traj = mdtraj.load(model_filepaths)

    def _mk_traj_alt(self):
        traj = mdtraj.load_pdb(self.model_filepaths[0])
        for model_filepath in self.model_filepaths[1:]:
            traj += mdtraj.load_pdb(model_filepath)
        self.traj = traj

    def rmsd(self):
        has_model_indices = self.df[self.df.has_model == True].index
        ca_atoms = [a.index for a in self.traj.topology.atoms if a.name == 'CA']
        rmsds = mdtraj.rmsd(self.traj, self.traj[0], atom_indices=ca_atoms, parallel=False)
        template_rmsds = [None] * len(self.templateids)
        for m,t in enumerate(has_model_indices):
            template_rmsds[t] = rmsds[m]

        self.df['rmsd'] = template_rmsds

    def rmsd_dist(self):
        warnings.warn('Not yet implemented.')
        pass

    def to_pickle(self, ofilepath):
        self.df.to_pickle(ofilepath)