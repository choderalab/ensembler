import os
import tempfile
import shutil
import mdtraj
from subprocess import Popen, PIPE
from ensembler.core import get_most_advanced_ensembler_modeling_stage, default_project_dirnames
from ensembler.core import model_filenames_by_ensembler_stage, get_valid_model_ids


class MolProbityValidation(object):
    def __init__(self,
                 targetid=None,
                 run_main=True
                 ):
        self.targetid = targetid
        self.models_target_dir = os.path.join(default_project_dirnames.models, self.targetid)
        self.molprobity_results_filepath = os.path.join(self.models_target_dir, 'molprobity-results')
        if run_main:
            self.tmp_models_dir = tempfile.mkdtemp()
            print self.tmp_models_dir
            try:
                self.ensembler_stage = get_most_advanced_ensembler_modeling_stage(self.targetid)
                self.model_structure_filename = model_filenames_by_ensembler_stage[self.ensembler_stage]
                self.valid_model_ids = get_valid_model_ids(
                    self.ensembler_stage,
                    self.targetid
                )
                self.copy_models_to_tmp_dir()
                self.run_molprobity()
                self.save_molprobity_results()
            finally:
                shutil.rmtree(self.tmp_models_dir)

    def get_models(self):
        valid_model_filepaths = []
        for filename in os.listdir(self.models_target_dir):
            if os.path.isdir(filename):
                model_filepath = os.path.join(self.models_target_dir, filename, self.model_structure_filename)
                if os.path.exists(model_filepath):
                    valid_model_filepaths.append(model_filepath)

    def copy_models_to_tmp_dir(self):
        """
        Also decompresses models and strips out waters.
        """
        for model_id in self.valid_model_ids:
            source_path = os.path.join(
                default_project_dirnames.models,
                self.targetid,
                model_id,
                self.model_structure_filename
            )

            dest_path = os.path.join(
                self.tmp_models_dir,
                model_id + '.pdb'
            )

            source_model_traj = mdtraj.load_pdb(source_path)
            protein_only_traj = source_model_traj.atom_slice(
                source_model_traj.top.select('protein')
            )
            protein_only_traj.save_pdb(dest_path)

    def run_molprobity(self):
        stdout, stderr = run_molprobity_oneline_analysis(self.tmp_models_dir)
        output_text = '\n'.join([stdout, stderr])
        self.molprobity_results = parse_molprobity_oneline_analysis_output(output_text)

    def save_molprobity_results(self):
        molprobity_results_tuples = [item for item in self.molprobity_results.iteritems()]
        molprobity_results_sorted = sorted(molprobity_results_tuples, key=lambda x: x[1])
        molprobity_results_str = '\n'.join(['{} {}'.format(*item) for item in molprobity_results_sorted])
        with open(self.molprobity_results_filepath, 'w') as molprobity_results_file:
            molprobity_results_file.write(molprobity_results_str + '\n')


def run_molprobity_oneline_analysis(dir_path):
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

    results = {}
    for results_line in results_lines:
        results_line_split = results_line.split(':')
        filename = results_line_split[0]
        targetid = filename[: filename.find('.pdb')]
        score_str = results_line_split[-2]
        try:
            score = float(score_str)
        except (TypeError, ValueError):
            continue
        results[targetid] = score
    return results
