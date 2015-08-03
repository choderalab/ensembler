import os
import Bio.SeqUtils
import mdtraj
from ensembler.core import check_project_toplevel_dir, default_project_dirnames
from ensembler.utils import set_loglevel
from ensembler.uniprot import get_uniprot_xml


class RenumberResidues(object):
    """Renumber residues according to the UniProt canonical sequence

    Requires the target id to be of the format '[UniProt mnemonic]_D[domain index]'

    Note that this will not work with targets with sequences which do not match the UniProt seq
    (e.g. those with mutations, insertions, deletions)
    """
    def __init__(self, targetid, project_dir='.', log_level=None):
        check_project_toplevel_dir()
        set_loglevel(log_level)
        self.targetid = targetid
        self.models_target_dir = os.path.join(default_project_dirnames.models, self.targetid)
        if not os.path.exists(self.models_target_dir):
            raise Exception('Model "{}" not found'.format(self.targetid))
        self.project_dir = project_dir
        self.uniprot_mnemonic = '_'.join(self.targetid.split('_')[0:2])
        self._get_models()
        self._get_model_seqs()
        self._get_uniprot_seq()
        self._find_seq_starts_and_ends()
        self._renumber_models()
        self._output_models()
        self._finish()

    def _get_models(self):
        self.model = {}
        root, dirnames, filenames = next(os.walk(self.models_target_dir))
        for dirname in dirnames:
            if 'implicit' in self.model and 'explicit' in self.model:
                break
            if 'implicit' not in self.model:
                implicit_model_filename = os.path.join(self.models_target_dir, dirname, 'implicit-refined.pdb.gz')
                if os.path.exists(implicit_model_filename):
                    self.model['implicit'] = mdtraj.load_pdb(implicit_model_filename)

            if 'explicit' not in self.model:
                explicit_model_filename = os.path.join(self.models_target_dir, dirname, 'explicit-refined.pdb.gz')
                if os.path.exists(explicit_model_filename):
                    self.model['explicit'] = mdtraj.load_pdb(explicit_model_filename)

    def _get_model_seqs(self):
        self.model_seq = ''.join([Bio.SeqUtils.seq1(r.name) for r in self.model['implicit'].top.residues])

    def _get_uniprot_seq(self):
        uniprot_query_string = 'mnemonic: {0}'.format(self.uniprot_mnemonic)
        self._uniprot_xml = get_uniprot_xml(uniprot_query_string)
        seq_rawtext = self._uniprot_xml.find('entry/sequence').text
        self.uniprot_seq = ''.join(seq_rawtext.split())

    def _find_seq_starts_and_ends(self):
        # 1-based inclusive numbering
        self.model_seq_start = self.uniprot_seq.index(self.model_seq) + 1
        self.model_seq_end = len(self.uniprot_seq) - self.uniprot_seq[::-1].index(self.model_seq[::-1])

    def _renumber_models(self):
        new_numbers = range(self.model_seq_start, self.model_seq_end+1)
        for key in self.model:
            for r, residue in enumerate(self.model[key].top.residues):
                if r < len(new_numbers):
                    residue.resSeq = new_numbers[r]
                else:
                    # ions and waters in explicit file should already be given chainID 'B' and numbered from 1
                    pass

    def _output_models(self):
        for key in self.model:
            ofilepath = os.path.join(self.models_target_dir, 'topol-renumbered-{0}.pdb'.format(key))
            self.model[key].save_pdb(ofilepath)

    def _finish(self):
        print('Done.')
