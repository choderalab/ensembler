from datetime import datetime
import ensembler
import os
import re
import pandas as pd
import yaml
import mdtraj
import gzip
from ensembler.core import logger
import warnings


class ProjectCounts(object):
    def __init__(self, targetid, project_dir='.', ofile_basepath='.', log_level=None):
        ensembler.utils.loglevel_setter(logger, log_level)
        self.targetid = targetid
        self.model_dir = os.path.join(ensembler.core.default_project_dirnames.models, self.targetid)
        self.project_dir = project_dir
        self.ofile_basepath = ofile_basepath
        self.df = pd.DataFrame()
        self._count_templates()
        self._count_models()
        self._count_uniques()
        self._count_implicit_refined()
        self._get_sequence_identities()

    def _count_templates(self):
        templateid = []
        root, dirnames, filenames = next(os.walk(self.model_dir))
        for dirname in dirnames:
            if re.match(ensembler.core.template_id_regex, dirname):
                templateid.append(dirname)
        self.df['templateid'] = templateid

    def _count_models(self):
        has_model = []
        for templateid in self.df.templateid:
            model_path = os.path.join(self.model_dir, templateid, 'model.pdb.gz')
            if os.path.exists(model_path):
                has_model.append(True)
            else:
                has_model.append(False)
        self.df['has_model'] = has_model

    def _count_uniques(self):
        unique = []
        for templateid in self.df.templateid:
            model_path = os.path.join(self.model_dir, templateid, 'unique_by_clustering')
            if os.path.exists(model_path):
                unique.append(True)
            else:
                unique.append(False)
        self.df['unique'] = unique

    def _count_implicit_refined(self):
        has_implicit_refined = []
        for templateid in self.df.templateid:
            model_path = os.path.join(self.model_dir, templateid, 'implicit-refined.pdb.gz')
            if os.path.exists(model_path):
                has_implicit_refined.append(True)
            else:
                has_implicit_refined.append(False)
        self.df['has_implicit_refined'] = has_implicit_refined

    def _get_sequence_identities(self):
        sequence_identity = []
        for templateid in self.df.templateid:
            aln_path = os.path.join(self.model_dir, templateid, 'alignment.pir')
            if os.path.exists(aln_path):
                with open(aln_path) as aln_file:
                    aln_text = aln_file.read().splitlines()
                aln = [aln_text[3], aln_text[6]]
                seqid = 100.0 * (sum([aa1 == aa2 for aa1, aa2 in zip(*aln)]) / float(len(aln[0])))
                sequence_identity.append(seqid)
            else:
                sequence_identity.append(None)
        self.df['sequence_identity'] = sequence_identity

    def _get_sequence_identities_old(self):   # TODO delete
        sequence_identity = []
        for templateid in self.df.templateid:
            seqid_path = os.path.join(self.model_dir, templateid, 'sequence-identity.txt')
            if os.path.exists(seqid_path):
                seqid = float(open(seqid_path).read().strip())
                sequence_identity.append(seqid)
            else:
                sequence_identity.append(None)
        self.df['sequence_identity'] = sequence_identity

    def save_df(self, ofilepath=None):
        if ofilepath is None:
            ofilepath = 'counts.csv'
            self.df.to_csv(ofilepath)

    def write_counts(self, ofilepath=None, seqid_range=None):
        if ofilepath is None:
            if seqid_range is None:
                ofilepath = os.path.join(self.ofile_basepath, 'counts.txt')
            else:
                ofilepath = os.path.join(self.ofile_basepath, 'counts-seqid{:.0f}-{:.0f}.txt'.format(seqid_range[0], seqid_range[1]))

        if seqid_range is None:
            df_selected = self.df
        else:
            df_selected = self.df[self.df.sequence_identity >= seqid_range[0]][self.df.sequence_identity < seqid_range[1]]

        counts = pd.Series(
            [
                seqid_range,
                len(df_selected),
                df_selected.has_model.sum(),
                df_selected.unique.sum(),
                df_selected.has_implicit_refined.sum(),
            ],
            index=[
                'sequence_identity_range',
                'templates',
                'models',
                'unique_models',
                'implicit_refined',
            ]
        )

        with open(ofilepath, 'w') as ofile:
            ofile.write(counts.to_string()+'\n')

    def write_attrition_rates(self, ofilepath=None, seqid_range=None):
        if ofilepath is None:
            if seqid_range is None:
                ofilepath = os.path.join(self.ofile_basepath, 'attrition_rates.txt')
            else:
                ofilepath = os.path.join(self.ofile_basepath, 'attrition_rates-seqid{:.0f}-{:.0f}.txt'.format(seqid_range[0], seqid_range[1]))

        if seqid_range is None:
            df_selected = self.df
        else:
            df_selected = self.df[self.df.sequence_identity >= seqid_range[0]][self.df.sequence_identity < seqid_range[1]]

        counts = [
            df_selected.has_model.sum(),
            df_selected.unique.sum(),
            df_selected.has_implicit_refined.sum(),
        ]

        rates = pd.Series(
            [
                seqid_range,
                1.0,
                counts[0] / float(len(df_selected)),
                counts[1] / float(counts[0]),
                counts[2] / float(counts[1]),
            ],
            index=[
                'sequence_identity_range',
                'templates',
                'models',
                'unique_models',
                'implicit_refined',
            ]
        )

        with open(ofilepath, 'w') as ofile:
            ofile.write(rates.to_string()+'\n')


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
        self.model_dir = os.path.join(ensembler.core.default_project_dirnames.models, self.targetid)
        self.project_dir = project_dir
        self._get_templateids_and_model_filepaths()
        self._get_unique_models()
        self._get_seqids()
        self._store_highest_seqid_model()
        self._mk_traj()

    def _get_templateids_and_model_filepaths(self):

        root, dirnames, filenames = os.walk(self.model_dir).next()

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

    def _get_seqids(self):
        seqid_filepath = os.path.join(self.model_dir, 'sequence-identities.txt')
        with open(seqid_filepath) as seqid_file:
            seqid_data = zip(*[line.split() for line in seqid_file.read().splitlines()])
        seqid_df = pd.DataFrame({
            'templateid': seqid_data[0],
            'seqid': seqid_data[1],
        })
        self.df = pd.merge(self.df, seqid_df, on='templateid')

    def _store_highest_seqid_model(self):
        models_sorted = self.df.sort('seqid').templateid
        for modelid in models_sorted:
            model_filepath = os.path.join(self.model_dir, modelid, 'model.pdb.gz')
            if os.path.exists(model_filepath):
                self.ref_modelid = modelid
                self.ref_model_filepath = model_filepath
                break
        self.ref_model_traj = mdtraj.load_pdb(self.ref_model_filepath)

    def rmsd(self):
        has_model_indices = self.df[self.df.has_model == True].index
        ca_atoms = [a.index for a in self.traj.topology.atoms if a.name == 'CA']
        rmsds = mdtraj.rmsd(self.traj, self.ref_model_traj, atom_indices=ca_atoms, parallel=False)
        template_rmsds = [None] * len(self.templateids)
        for m,t in enumerate(has_model_indices):
            template_rmsds[t] = rmsds[m]

        self.df['rmsd'] = template_rmsds

    def rmsd_dist(self):
        warnings.warn('Not yet implemented.')
        pass

    def to_csv(self, ofilepath=None):
        if ofilepath is None:
            ofilepath = os.path.join(self.model_dir, 'rmsds.csv')
        self.df.to_csv(ofilepath)