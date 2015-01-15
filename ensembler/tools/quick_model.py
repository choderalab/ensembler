import os
import pandas as pd
import ensembler
import ensembler.initproject
import ensembler.modeling
import ensembler.refinement
import ensembler.packaging


class QuickParentClass(object):
    pass


class QuickModel(QuickParentClass):
    def __init__(self, targetids, templateids=None, template_seqid_cutoff=None, loopmodel=False, package_for_fah=False, nfahclones=None):
        """
        Run this after having set up targets and templates with the appropriate ensembler commands.

        :param targetids: list of str
        :param templateids: list of list of str (one list of templateids for each target)
        :param template_seqid_cutoff:
        :param package_for_fah:
        :param nfahclones:
        """
        self.targetids = targetids
        self.loopmodel = loopmodel
        if templateids:
            self.templateids = templateids
        elif template_seqid_cutoff:
            self.align_all_templates()
            self.templateids = self.select_templates_based_on_seqid_cutoff(seqid_cutoff=template_seqid_cutoff)
        else:
            self.align_all_templates()
            self.templateids = self.select_templates_based_on_seqid_cutoff()

        print self.templateids # TODO
        self.model(package_for_fah=package_for_fah, nfahclones=nfahclones) # TODO

    def model(self, package_for_fah=False, nfahclones=None):
        # ensembler.initproject.GatherTargetsFromUniProt()
        for t in range(len(self.targetids)):
            if self.loopmodel:
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

    def align_all_templates(self):
        ensembler.modeling.align_targets_and_templates(process_only_these_targets=self.targetids)

    def select_templates_based_on_seqid_cutoff(self, seqid_cutoff=None):
        """
        If seqid_cutoff is not specified, it will be chosen interactively.
        :param seqid_cutoff:
        :return:
        """
        for targetid in self.targetids:
            seqid_filepath = os.path.join(ensembler.core.default_project_dirnames.models, targetid, 'sequence-identities.txt')
            with open(seqid_filepath) as seqid_file:
                seqids_list = [line.split() for line in seqid_file.readlines()]
                templateids = [i[0] for i in seqids_list]
                seqids = [float(i[1]) for i in seqids_list]
            df = pd.DataFrame({
                'templateids': templateids,
                'seqids': seqids,
                })
            if seqid_cutoff:
                selected_templates = df[df.seqids > seqid_cutoff]
                templateids = list(selected_templates.templateids)
            else:
                for cutoff in range(0, 101, 10):
                    ntemplates = len(df[df.seqids > cutoff])
                    print 'Number of templates with seqid > %d: %d' % (cutoff, ntemplates)
                seqid_cutoff_chosen = False

                while not seqid_cutoff_chosen:
                    seqid_cutoff = float(raw_input('Choose a sequence identity cutoff (confirm at next step): '))
                    ntemplates = len(df[df.seqids > seqid_cutoff])
                    print 'Number of templates chosen by sequence identity cutoff (%.1f): %d' % (seqid_cutoff, ntemplates)
                    confirm_seqid_cutoff = raw_input('Use this sequence identity cutoff? (y|N): ')
                    if confirm_seqid_cutoff.lower() in ['y', 'yes']:
                        seqid_cutoff_chosen = True

                selected_templates = df[df.seqids > seqid_cutoff]
                templateids = list(selected_templates.templateids)

        return templateids


class QuickEnsemble(QuickParentClass):
    def __init__(self, input_ids, uniprot_entry_name, pdbids=None, chainids=None, template_seqid_cutoff=None, package_for_fah=False, nfahclones=None):
        """
        Includes target and template retrieval.

        :param input_ids: {uniprot_ac: {pdbid: [chainid, chainid]}, uniprot_ac: {}}
        :param uniprot_entry_name: str
        :param pdbids: list of str [pdbid, pdbid]
        :param chainids: list of list of str [[chainid, chainid], []]
        :param template_seqid_cutoff:
        :param package_for_fah:
        :param nfahclones:
        """
        pass