import ensembler
import ensembler.initproject
import ensembler.modeling


class QuickModel:
    def __init__(self, targetids, templateids=None, seqid_cutoff=None):
        self.targetids = targetids
        self.templateids = templateids
        self.seqid_cutoff = seqid_cutoff
        self.main()

    def main(self):
        # ensembler.initproject.GatherTargetsFromUniProt()
        ensembler.modeling.align_targets_and_templates(targets=self.targetids, templates=self.templateids)
        ensembler.modeling.extract_template_structures_from_pdb_files(templates=self.templateids)
        ensembler.modeling.loopmodel_templates(templates=self.templateids)