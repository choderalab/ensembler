import ensembler
import ensembler.initproject
import ensembler.modeling


class QuickModel:
    def __init__(self, targetids, templateids=None, template_seqid_cutoff=None):
        self.targetids = targetids
        if templateids:
            self.templateids = templateids
        elif template_seqid_cutoff:
            self.template_seqid_cutoff = template_seqid_cutoff
            templateids = self.select_templates_based_on_seqid_cutoff(template_seqid_cutoff)
        self.main()

    def main(self):
        # ensembler.initproject.GatherTargetsFromUniProt()
        ensembler.modeling.model_template_loops(targets=self.targetids, templates=self.templateids)
        ensembler.modeling.align_targets_and_templates(targets=self.targetids, templates=self.templateids)

    def select_templates_based_on_seqid_cutoff(self, seqid_cutoff):
        raise NotImplementedError
        return templateids
