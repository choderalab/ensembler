import ensembler
import ensembler.initproject
import ensembler.modeling
import ensembler.refinement
import ensembler.packaging


class QuickModel:
    def __init__(self, targetids, templateids=None, template_seqid_cutoff=None, package_for_fah=False, nfahclones=None):
        """
        :param targetids: list of str
        :param templateids: list of list of str (one list of templateids for each target)
        :param template_seqid_cutoff:
        :param package_for_fah:
        :param nfahclones:
        :return:
        """
        self.targetids = targetids
        if templateids:
            self.templateids = templateids
        elif template_seqid_cutoff:
            self.align_without_specified_templates()
            self.templateids = self.select_templates_based_on_seqid_cutoff(template_seqid_cutoff=template_seqid_cutoff)
        else:
            self.align_without_specified_templates()
            self.templateids = self.select_templates_based_on_seqid_cutoff()

        print self.templateids
        # self.model(package_for_fah=package_for_fah, nfahclones=nfahclones)

    def model(self, package_for_fah=False, nfahclones=None):
        # ensembler.initproject.GatherTargetsFromUniProt()
        for t in range(len(self.targetids)):
            ensembler.modeling.model_template_loops(process_only_these_templates=self.templateids[t])
            ensembler.modeling.align_targets_and_templates(process_only_these_targets=self.targetids, process_only_these_templates=self.templateids[t])
            ensembler.modeling.build_models(process_only_these_targets=self.targetids, process_only_these_templates=self.templateids[t])
            ensembler.modeling.cluster_models(process_only_these_targets=self.targetids)
            ensembler.refinement.refine_implicit_md(process_only_these_targets=self.targetids, process_only_these_templates=self.templateids[t])
            ensembler.refinement.solvate_models(process_only_these_targets=self.targetids, process_only_these_templates=self.templateids[t])
            ensembler.refinement.determine_nwaters(process_only_these_targets=self.targetids, process_only_these_templates=self.templateids[t])
            ensembler.refinement.refine_explicitMD(process_only_these_targets=self.targetids, process_only_these_templates=self.templateids[t])
            if package_for_fah:
                if not nfahclones:
                    nfahclones = 1
                ensembler.packaging.package_for_fah(process_only_these_targets=self.targetids, nclones=nfahclones)

    def align_without_specified_templates(self):
        ensembler.modeling.align_targets_and_templates(process_only_these_targets=self.targetids)

    def select_templates_based_on_seqid_cutoff(self, template_seqid_cutoff=None):
        if template_seqid_cutoff:
            pass
            # use specified seqid_cutoff to select templateids
        else:
            pass
            # select interactively

        raise NotImplementedError
        return templateids
