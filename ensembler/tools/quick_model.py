import os
import warnings
import pandas as pd
import ensembler
import ensembler.initproject
import ensembler.modeling
import ensembler.refinement
import ensembler.packaging


class QuickParentClass(object):
    def __init__(self, targetid=None, templateids=None):
        self.targetid = targetid
        self.templateids = templateids

    def _model(self, loopmodel=True, package_for_fah=False, nfahclones=None):
        if loopmodel:
            ensembler.modeling.model_template_loops(process_only_these_templates=self.templateids)
        ensembler.modeling.align_targets_and_templates(process_only_these_targets=[self.targetid], process_only_these_templates=self.templateids)
        ensembler.modeling.build_models(process_only_these_targets=[self.targetid], process_only_these_templates=self.templateids)
        ensembler.modeling.cluster_models(process_only_these_targets=[self.targetid])
        ensembler.refinement.refine_implicit_md(process_only_these_targets=[self.targetid], process_only_these_templates=self.templateids)
        ensembler.refinement.solvate_models(process_only_these_targets=[self.targetid], process_only_these_templates=self.templateids)
        ensembler.refinement.determine_nwaters(process_only_these_targets=[self.targetid], process_only_these_templates=self.templateids)
        ensembler.refinement.refine_explicitMD(process_only_these_targets=[self.targetid], process_only_these_templates=self.templateids)
        if package_for_fah:
            if not nfahclones:
                nfahclones = 1
            ensembler.packaging.package_for_fah(process_only_these_targets=[self.targetid], nclones=nfahclones)

    def _align_all_templates(self):
        ensembler.modeling.align_targets_and_templates(process_only_these_targets=[self.targetid])

    def _select_templates_based_on_seqid_cutoff(self, seqid_cutoff=None):
        """
        If seqid_cutoff is not specified, it will be chosen interactively.
        :param seqid_cutoff:
        :return:
        """
        seqid_filepath = os.path.join(ensembler.core.default_project_dirnames.models, self.targetid, 'sequence-identities.txt')
        with open(seqid_filepath) as seqid_file:
            seqids_list = [line.split() for line in seqid_file.readlines()]
            templateids = [i[0] for i in seqids_list]
            seqids = [float(i[1]) for i in seqids_list]

        df = pd.DataFrame({
            'templateids': templateids,
            'seqids': seqids,
            })

        if not seqid_cutoff:
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


class QuickModel(QuickParentClass):
    def __init__(self, targetid, templateids=None, template_seqid_cutoff=None, loopmodel=True, package_for_fah=False, nfahclones=None):
        """
        Run this after having set up targets and templates with the appropriate ensembler commands.

        :param targetid: str
        :param templateids: list of str
        :param template_seqid_cutoff:
        :param package_for_fah:
        :param nfahclones:
        """
        super(QuickModel, self).__init__()
        self.targetid = targetid

        all_templates_resolved_seq, all_templates_full_seq = ensembler.core.get_templates()
        if len(all_templates_resolved_seq) == 0:
            warnings.warn('No templates found. Exiting.')
            return

        if templateids:
            self.templateids = templateids
        else:
            self._align_all_templates()
            self.templateids = self._select_templates_based_on_seqid_cutoff(seqid_cutoff=template_seqid_cutoff)

        if len(self.templateids) == 0:
            warnings.warn('No templates selected. Exiting.')
            return

        self._model(loopmodel=loopmodel, package_for_fah=package_for_fah, nfahclones=nfahclones)


class QuickEnsemble(QuickParentClass):
    def __init__(self, target_uniprot_entry_name, uniprot_domain_regex, pdbids=None, chainids=None, template_seqid_cutoff=None, loopmodel=False, package_for_fah=False, nfahclones=None, project_dir=None, structure_dirs=None):
        """
        Includes target and template retrieval.

        :param target_uniprot_entry_name: str
        :param uniprot_domain_regex: str
        :param pdbids: list of str [pdbid, pdbid]
        :param chainids: list of list of str [[chainid, chainid], []]
        :param template_seqid_cutoff:
        :param package_for_fah:
        :param nfahclones:
        """
        super(QuickEnsemble, self).__init__()
        if project_dir:
            os.chdir(project_dir)
        ensembler.initproject.InitProject('.')

        uniprot_query_string = 'mnemonic: %s' % target_uniprot_entry_name
        gather_targets_obj = ensembler.initproject.GatherTargetsFromUniProt(uniprot_query_string, uniprot_domain_regex=uniprot_domain_regex)
        self.targetid = [target.id for target in gather_targets_obj.targets]

        if pdbids:
            ensembler.initproject.gather_templates_from_pdb(pdbids, uniprot_domain_regex, chainids=chainids, structure_dirs=structure_dirs)
            templates_resolved_seq, templates_full_seq = ensembler.core.get_templates()
            self.templateids = [t.id for t in templates_resolved_seq]
        else:
            template_uniprot_query_string = 'domain: "%s" AND reviewed: yes' % uniprot_domain_regex
            ensembler.initproject.gather_templates_from_uniprot(template_uniprot_query_string, uniprot_domain_regex, structure_dirs=structure_dirs)
            self._align_all_templates()
            self.templateids = self._select_templates_based_on_seqid_cutoff(seqid_cutoff=template_seqid_cutoff)

        self._model(loopmodel=loopmodel, package_for_fah=package_for_fah, nfahclones=nfahclones)