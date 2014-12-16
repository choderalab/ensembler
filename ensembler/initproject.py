import gzip
import json
import sys
import os
import shutil
import subprocess
import traceback
import tempfile
import datetime
from lxml import etree
import yaml
import Bio.SeqUtils
import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import ensembler
import ensembler.TargetExplorer
import ensembler.UniProt
import ensembler.PDB
import ensembler.version
from ensembler.utils import file_exists_and_not_empty
from ensembler.core import mpistate, logger


class TemplateData:
    def __init__(self, pdbid=None, chainid=None, templateid=None, resolved_seq=None, resolved_pdbresnums=None, full_seq=None, full_pdbresnums=None):
        self.pdbid = pdbid
        self.chainid = chainid
        self.templateid = templateid
        self.resolved_seq = resolved_seq
        self.resolved_pdbresnums = resolved_pdbresnums
        self.full_seq = full_seq
        self.full_pdbresnums = full_pdbresnums


class LoopmodelOutput:
    def __init__(self, output_text=None, loopmodel_exception=None, exception=None, trbk=None, successful=False, no_missing_residues=False):
        self.output_text = output_text
        self.exception = exception
        self.loopmodel_exception = loopmodel_exception
        self.traceback = trbk
        self.successful = successful
        self.no_missing_residues = no_missing_residues


class InitProject(object):
    def __init__(self, project_toplevel_dir, run_main=True):
        """Initialize Ensembler project within the given directory. Creates
        necessary subdirectories and a project metadata .yaml file.
        :param project_toplevel_dir: str
        """
        self.project_toplevel_dir = project_toplevel_dir
        if run_main:
            self._init_project()

    @ensembler.utils.notify_when_done
    def _init_project(self):
        self._create_project_dirs()
        self._write_init_metadata()

    def _create_project_dirs(self):
        ensembler.utils.create_dir(os.path.join(self.project_toplevel_dir, ensembler.core.default_project_dirnames.targets))
        ensembler.utils.create_dir(os.path.join(self.project_toplevel_dir, ensembler.core.default_project_dirnames.templates))
        ensembler.utils.create_dir(os.path.join(self.project_toplevel_dir, ensembler.core.default_project_dirnames.structures))
        ensembler.utils.create_dir(os.path.join(self.project_toplevel_dir, ensembler.core.default_project_dirnames.models))
        ensembler.utils.create_dir(os.path.join(self.project_toplevel_dir, ensembler.core.default_project_dirnames.packaged_models))
        ensembler.utils.create_dir(os.path.join(self.project_toplevel_dir, ensembler.core.default_project_dirnames.structures_pdb))
        ensembler.utils.create_dir(os.path.join(self.project_toplevel_dir, ensembler.core.default_project_dirnames.structures_sifts))
        ensembler.utils.create_dir(os.path.join(self.project_toplevel_dir, ensembler.core.default_project_dirnames.templates_structures_resolved))
        ensembler.utils.create_dir(os.path.join(self.project_toplevel_dir, ensembler.core.default_project_dirnames.templates_structures_modeled_loops))

    def _write_init_metadata(self):
        project_metadata = ensembler.core.ProjectMetadata(project_stage='init', project_toplevel_dir=self.project_toplevel_dir)
        init_metadata = self._gen_init_metadata()
        project_metadata.add_data(init_metadata)
        project_metadata.write()

    def _gen_init_metadata(self):
        datestamp = ensembler.core.get_utcnow_formatted()
        metadata_dict = {
            'datestamp': datestamp,
            'init_path': os.path.abspath(self.project_toplevel_dir),
            'python_version': sys.version.split('|')[0].strip(),
            'python_full_version': ensembler.core.literal_str(sys.version),
            'ensembler_version': ensembler.version.short_version,
            'ensembler_commit': ensembler.version.git_revision
        }
        return metadata_dict


class GatherTargets(object):
    def __init__(self):
        self.manual_overrides = ensembler.core.ManualOverrides()

    def _gen_gather_targets_metadata(self, ntargets, additional_metadata=None):
        if additional_metadata is None:
            additional_metadata = {}
        datestamp = ensembler.core.get_utcnow_formatted()
        metadata = {
            'datestamp': datestamp,
            'ntargets': '%d' % ntargets,
            'python_version': sys.version.split('|')[0].strip(),
            'python_full_version': ensembler.core.literal_str(sys.version),
            'ensembler_version': ensembler.version.short_version,
            'ensembler_commit': ensembler.version.git_revision,
        }
        metadata.update(additional_metadata)
        return metadata


class GatherTargetsFromTargetExplorer(GatherTargets):
    def __init__(self, dbapi_uri, search_string='', run_main=True):
        super(GatherTargetsFromTargetExplorer, self).__init__()
        self.dbapi_uri = dbapi_uri
        self.search_string = search_string
        if run_main:
            self._gather_targets()

    @ensembler.utils.notify_when_done
    def _gather_targets(self):
        targetexplorer_json = ensembler.TargetExplorer.get_targetexplorer_json(self.dbapi_uri, self.search_string, return_data='domain_seqs,seqs')
        self.targets = self._extract_targets_from_json(targetexplorer_json)
        fasta_ofilepath = os.path.join(ensembler.core.default_project_dirnames.targets, 'targets.fa')
        Bio.SeqIO.write(self.targets, fasta_ofilepath, 'fasta')
        self._write_metadata()

    def _extract_targets_from_json(self, targetexplorer_json):
        """
        :param targetexplorer_json: dict
        :return: list of BioPython SeqRecord objects
        """
        targets = []
        for target in targetexplorer_json['results']:
            for target_domain in target['domains']:
                targetid = target_domain.get('targetid')
                targetseq = target_domain.get('sequence')
                # domain span override
                if targetid in self.manual_overrides.target.domain_spans:
                    fullseq = target.get('sequence')
                    start, end = [int(x) - 1 for x in self.manual_overrides.target.domain_spans[targetid].split('-')]
                    targetseq = fullseq[start:end + 1]
                target = SeqRecord(Seq(targetseq), id=targetid, description=targetid)
                targets.append(target)

        return targets

    def _write_metadata(self):
        targetexplorer_metadata = gen_targetexplorer_metadata(self.dbapi_uri, self.search_string)
        gather_targets_metadata = self._gen_gather_targets_metadata(len(self.targets), additional_metadata=targetexplorer_metadata)
        project_metadata = ensembler.core.ProjectMetadata(project_stage='gather_targets')
        project_metadata.add_data(gather_targets_metadata)
        project_metadata.write()


def gen_targetexplorer_metadata(dbapi_uri, search_string):
    db_metadata = ensembler.TargetExplorer.get_targetexplorer_metadata(dbapi_uri)
    metadata = {
        'method': 'TargetExplorer',
        'gather_from_targetexplorer': {
            'db_uniprot_query_string': str(db_metadata.get('uniprot_query_string')),
            'db_uniprot_domain_regex': str(db_metadata.get('uniprot_domain_regex')),
            'search_string': search_string,
            'dbapi_uri': dbapi_uri,
        }
    }
    return metadata


class GatherTargetsFromUniProt(GatherTargets):
    def __init__(self, uniprot_query_string, uniprot_domain_regex='', run_main=True):
        super(GatherTargetsFromUniProt, self).__init__()
        self.uniprot_query_string = uniprot_query_string
        self.uniprot_domain_regex = uniprot_domain_regex
        if run_main:
            self._gather_targets()

    @ensembler.utils.notify_when_done
    def _gather_targets(self):
        logger.info('Querying UniProt web server...')
        uniprotxml = ensembler.UniProt.get_uniprot_xml(self.uniprot_query_string)
        logger.info('Number of entries returned from initial UniProt search: %r\n' % len(uniprotxml))
        log_unique_domain_names(self.uniprot_query_string, uniprotxml)
        if self.uniprot_domain_regex is not None:
            log_unique_domain_names_selected_by_regex(self.uniprot_domain_regex, uniprotxml)
        fasta_ofilepath = os.path.join(ensembler.core.default_project_dirnames.targets, 'targets.fa')
        self.targets = self._extract_targets_from_uniprot_xml(uniprotxml)
        Bio.SeqIO.write(self.targets, fasta_ofilepath, 'fasta')
        self._write_metadata()

    def _extract_targets_from_uniprot_xml(self, uniprotxml):
        targets = []
        for entry in uniprotxml.findall('entry'):
            entry_name = entry.find('name').text
            fullseq = ensembler.core.sequnwrap(entry.find('sequence').text)
            if self.uniprot_domain_regex is not None:
                selected_domains = entry.xpath(
                    'feature[@type="domain"][match_regex(@description, "%s")]' % self.uniprot_domain_regex,
                    extensions={(None, 'match_regex'): ensembler.core.xpath_match_regex_case_sensitive}
                )
            else:
                selected_domains = entry.findall('feature[@type="domain"]')

            domain_iter = 0
            for domain in selected_domains:
                targetid = '%s_D%d' % (entry_name, domain_iter)
                # domain span override
                if targetid in self.manual_overrides.target.domain_spans:
                    start, end = [int(x) - 1 for x in self.manual_overrides.target.domain_spans[targetid].split('-')]
                else:
                    start, end = [int(domain.find('location/begin').get('position')) - 1,
                                  int(domain.find('location/end').get('position')) - 1]
                targetseq = fullseq[start:end + 1]
                targets.append(SeqRecord(Seq(targetseq), id=targetid, description=targetid))
                domain_iter += 1

        return targets

    def _write_metadata(self):
        uniprot_metadata = gen_uniprot_metadata(self.uniprot_query_string, self.uniprot_domain_regex)
        gather_targets_metadata = self._gen_gather_targets_metadata(len(self.targets), additional_metadata=uniprot_metadata)
        project_metadata = ensembler.core.ProjectMetadata(project_stage='gather_targets')
        project_metadata.add_data(gather_targets_metadata)
        project_metadata.write()


def gen_uniprot_metadata(uniprot_query_string, uniprot_domain_regex):
    metadata = {
        'method': 'UniProt',
        'gather_from_uniprot': {
            'uniprot_query_string': uniprot_query_string,
            'uniprot_domain_regex': uniprot_domain_regex,
        }
    }
    return metadata


def log_unique_domain_names(uniprot_query_string, uniprotxml):
    if 'domain:' in uniprot_query_string:
        # First extract the domain selection
        # Example query string: 'domain:"Protein kinase" AND reviewed:yes'
        # Can assume that the domain selection will be bounded by double-quotes
        query_string_split = uniprot_query_string.split('"')
        query_string_domain_selection = query_string_split[query_string_split.index('domain:') + 1]
        uniprot_query_string_domains = uniprotxml.xpath(
            'entry/feature[@type="domain"][match_regex(@description, "%s")]' % query_string_domain_selection,
            extensions={
                (None, 'match_regex'): ensembler.core.xpath_match_regex_case_insensitive
            }
        )
        uniprot_unique_domain_names = set([domain.get('description') for domain in uniprot_query_string_domains])
        logger.info('Set of unique domain names selected by the domain selector \'%s\' during the initial UniProt search:\n%s\n'
              % (query_string_domain_selection, uniprot_unique_domain_names))

    else:
        uniprot_domains = uniprotxml.xpath('entry/feature[@type="domain"]')
        uniprot_unique_domain_names = set([domain.get('description') for domain in uniprot_domains])
        logger.info('Set of unique domain names returned from the initial UniProt search using the query string \'%s\':\n%s\n'
              % (uniprot_query_string, uniprot_unique_domain_names))


def log_unique_domain_names_selected_by_regex(uniprot_domain_regex, uniprotxml):
    regex_matched_domains = uniprotxml.xpath(
        'entry/feature[@type="domain"][match_regex(@description, "%s")]' % uniprot_domain_regex,
        extensions={(None, 'match_regex'): ensembler.core.xpath_match_regex_case_sensitive}
    )
    regex_matched_domains_unique_names = set([domain.get('description') for domain in regex_matched_domains])
    logger.info('Unique domain names selected after searching with the case-sensitive regex string \'%s\':\n%s\n'
        % (uniprot_domain_regex, regex_matched_domains_unique_names))


@ensembler.utils.notify_when_done
def gather_templates_from_targetexplorer(dbapi_uri, search_string='', structure_dirs=None, loopmodel=True, overwrite_structures=False):
    """Gather protein template data from a TargetExplorer DB network API.
    Pass the URI for the database API and a search string.
    The search string uses SQLAlchemy syntax and standard TargetExplorer
    frontend data fields.
    Example:
    dbapi_uri='http://plfah2.mskcc.org/kinomeDBAPI'
    search_string='species="Human"'

    To select all domains within the database:
    search_string=''
    """
    manual_overrides = ensembler.core.ManualOverrides()
    templates_json = get_targetexplorer_templates_json(dbapi_uri, search_string)
    selected_pdbchains = extract_template_pdbchains_from_targetexplorer_json(templates_json, manual_overrides=manual_overrides)
    for pdbchain in selected_pdbchains:
        get_pdb_and_sifts_files(pdbchain['pdbid'], structure_dirs)

    selected_templates = extract_template_pdb_chain_residues(selected_pdbchains)
    write_template_seqs_to_fasta_file(selected_templates)
    extract_template_structures_from_pdb_files(selected_templates, overwrite_structures=overwrite_structures)
    if loopmodel:
        missing_residues = pdbfix_templates(selected_templates, overwrite_structures=overwrite_structures)
        loopmodel_templates(selected_templates, missing_residues, overwrite_structures=overwrite_structures)
    write_gather_templates_from_targetexplorer_metadata(search_string, dbapi_uri, len(selected_templates), structure_dirs, loopmodel)


@ensembler.utils.notify_when_done
def gather_templates_from_uniprot(uniprot_query_string, uniprot_domain_regex, structure_dirs=None, loopmodel=True, overwrite_structures=False):
    """# Searches UniProt for a set of template proteins with a user-defined
    query string, then saves IDs, sequences and structures."""
    manual_overrides = ensembler.core.ManualOverrides()
    selected_pdbchains = None
    if mpistate.rank == 0:
        uniprotxml = ensembler.UniProt.get_uniprot_xml(uniprot_query_string)
        log_unique_domain_names(uniprot_query_string, uniprotxml)
        if uniprot_domain_regex is not None:
            log_unique_domain_names_selected_by_regex(uniprot_domain_regex, uniprotxml)

        selected_pdbchains = extract_template_pdbchains_from_uniprot_xml(uniprotxml, uniprot_domain_regex, manual_overrides)
        for pdbchain in selected_pdbchains:
            get_pdb_and_sifts_files(pdbchain['pdbid'], structure_dirs)

    selected_pdbchains = mpistate.comm.bcast(selected_pdbchains, root=0)

    selected_templates = extract_template_pdb_chain_residues(selected_pdbchains)
    write_template_seqs_to_fasta_file(selected_templates)
    extract_template_structures_from_pdb_files(selected_templates, overwrite_structures=overwrite_structures)
    if loopmodel:
        missing_residues = pdbfix_templates(selected_templates, overwrite_structures=overwrite_structures)
        print 'Starting loopmodeling'
        loopmodel_templates(selected_templates, missing_residues, overwrite_structures=overwrite_structures)
    write_gather_templates_from_uniprot_metadata(uniprot_query_string, uniprot_domain_regex, len(selected_templates), structure_dirs, loopmodel)


def get_targetexplorer_templates_json(dbapi_uri, search_string):
    """
    :param dbapi_uri: str
    :param search_string: str
    :return: list containing nested lists and dicts
    """
    targetexplorer_json = None
    if mpistate.rank == 0:
        targetexplorer_jsonstr = ensembler.TargetExplorer.query_targetexplorer(
            dbapi_uri, search_string, return_data='pdb_data'
        )
        targetexplorer_json = json.loads(targetexplorer_jsonstr)
    targetexplorer_json = mpistate.comm.bcast(targetexplorer_json, root=0)
    return targetexplorer_json


def extract_template_pdbchains_from_targetexplorer_json(targetexplorer_json, manual_overrides):
    selected_pdbchains = None
    if mpistate.rank == 0:
        selected_pdbchains = []
        for target in targetexplorer_json['results']:
            entry_name = target['entry_name']
            for pdb_data in target['pdbs']:
                pdbid = pdb_data['pdbid']
                for pdbchain_data in pdb_data['pdbchains']:
                    pdbchain_data['entry_name'] = entry_name
                    pdbchain_data['pdbid'] = pdbid
                    targetid = '%s_D%s' % (entry_name, pdbchain_data['domainid'])
                    templateid = '%s_%s_%s' % (targetid, pdbid, pdbchain_data['chainid'])
                    pdbchain_data['templateid'] = templateid
                    pdbchain_data['domain_span'] = [int(pdbchain_data['seq_begin']), int(pdbchain_data['seq_end'])]
                    # manual overrides
                    domain_len = pdbchain_data['seq_end'] - pdbchain_data['seq_begin'] + 1
                    if manual_overrides.template.min_domain_len is not None and domain_len < manual_overrides.template.min_domain_len:
                        continue
                    if manual_overrides.template.max_domain_len is not None and domain_len > manual_overrides.template.max_domain_len:
                        continue
                    if pdbid in manual_overrides.template.skip_pdbs:
                        continue
                    if targetid in manual_overrides.template.domain_spans:
                        pdbchain_data['domain_span'] = [int(x) for x in manual_overrides.template.domain_spans[targetid].split('-')]
                    selected_pdbchains.append(pdbchain_data)
    selected_pdbchains = mpistate.comm.bcast(selected_pdbchains, root=0)
    return selected_pdbchains


def attempt_symlink_structure_files(pdbid, project_structures_dir, structure_dirs, structure_type='pdb'):
    project_structure_filepath = os.path.join(project_structures_dir, structure_type, pdbid + structure_type_file_extension_mapper[structure_type])
    for structure_dir in structure_dirs:
        structure_filepath = os.path.join(structure_dir, pdbid + structure_type_file_extension_mapper[structure_type])
        if os.path.exists(structure_filepath):
            if file_exists_and_not_empty(structure_filepath) > 0:
                if os.path.exists(project_structure_filepath):
                    os.remove(project_structure_filepath)
                os.symlink(structure_filepath, project_structure_filepath)
                break


def download_structure_file(pdbid, project_structure_filepath, structure_type='pdb'):
    if structure_type == 'pdb':
        download_pdb_file(pdbid, project_structure_filepath)
    elif structure_type == 'sifts':
        download_sifts_file(pdbid, project_structure_filepath)


def download_pdb_file(pdbid, project_pdb_filepath):
    logger.info('Downloading PDB file for: %s' % pdbid)
    pdbgz_page = ensembler.PDB.retrieve_pdb(pdbid, compressed='yes')
    with open(project_pdb_filepath, 'w') as pdbgz_file:
        pdbgz_file.write(pdbgz_page)


def download_sifts_file(pdbid, project_sifts_filepath):
    logger.info('Downloading sifts file for: %s', pdbid)
    sifts_page = ensembler.PDB.retrieve_sifts(pdbid)
    with gzip.open(project_sifts_filepath, 'wb') as project_sifts_file:
        project_sifts_file.write(sifts_page)


def get_pdb_and_sifts_files(pdbid, structure_dirs=None):
    if type(structure_dirs) != list:
        structure_dirs = []
    project_structures_dir = 'structures'
    for structure_type in ['pdb', 'sifts']:
        project_structure_filepath = os.path.join(project_structures_dir, structure_type, pdbid + structure_type_file_extension_mapper[structure_type])
        if not file_exists_and_not_empty(project_structure_filepath):
            attempt_symlink_structure_files(pdbid, project_structures_dir, structure_dirs, structure_type=structure_type)
            if not os.path.exists(project_structure_filepath):
                download_structure_file(pdbid, project_structure_filepath, structure_type=structure_type)


def extract_template_pdb_chain_residues(selected_pdbchains):
    selected_templates = None
    if mpistate.rank == 0:
        logger.info('Extracting residues from PDB chains...')
        selected_templates = []
        for pdbchain in selected_pdbchains:
            extracted_pdb_template_seq_data = extract_pdb_template_seq(pdbchain)
            if extracted_pdb_template_seq_data is not None:
                selected_templates.append(extracted_pdb_template_seq_data)
        logger.info('%d templates selected.\n' % len(selected_templates))
    selected_templates = mpistate.comm.bcast(selected_templates, root=0)
    return selected_templates


def parse_sifts_xml(sifts_filepath):
    with gzip.open(sifts_filepath, 'rb') as sifts_file:
        parser = etree.XMLParser(huge_tree=True)
        siftsxml = etree.parse(sifts_file, parser).getroot()

    return siftsxml


def add_pdb_modified_xml_tags_to_residues(siftsxml):
    """
    Adds "PDB modified" tags to certain phosphorylated residue types, which sometimes do not have such tags in the SIFTS file.
    known cases: 4BCP, 4BCG, 4I5C, 4IVB, 4IAC
    The passed XML object is modified in-place.
    :param siftsxml:
    :return:
    """
    modified_residues = []
    modified_residues += siftsxml.findall('entity/segment/listResidue/residue[@dbResName="TPO"]')
    modified_residues += siftsxml.findall('entity/segment/listResidue/residue[@dbResName="PTR"]')
    modified_residues += siftsxml.findall('entity/segment/listResidue/residue[@dbResName="SEP"]')
    for mr in modified_residues:
        if mr is None:
            continue
        residue_detail_modified = etree.Element('residueDetail')
        residue_detail_modified.set('dbSource', 'MSD')
        residue_detail_modified.set('property', 'Annotation')
        residue_detail_modified.text = 'PDB\n          modified'
        mr.append(residue_detail_modified)


def extract_pdb_template_seq(pdbchain):
    """Extract data from PDB chain"""
    templateid = pdbchain['templateid']
    chainid = pdbchain['chainid']
    pdbid = pdbchain['pdbid']
    domain_span = pdbchain['domain_span'] # UniProt coords

    sifts_filepath = os.path.join('structures', 'sifts', pdbid + '.xml.gz')
    siftsxml = parse_sifts_xml(sifts_filepath)

    add_pdb_modified_xml_tags_to_residues(siftsxml)

    # domain_span_sifts_coords = [
    #     int(siftsxml.find('entity/segment/listResidue/residue/crossRefDb[@dbSource="UniProt"][@dbChainId="%s"][@dbResNum="%d"]/..' % (chainid, domain_span[0])).get('dbResNum')),
    #     int(siftsxml.find('entity/segment/listResidue/residue/crossRefDb[@dbSource="UniProt"][@dbChainId="%s"][@dbResNum="%d"]/..' % (chainid, domain_span[1])).get('dbResNum')),
    # ]

    # TODO once working, check whether this fills in all loops
    # An alternative approach would be to just take the UniProt sequence specified by the domain span
    selected_residues = siftsxml.xpath(
        'entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"]'
        '[../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]]' % (chainid, domain_span[0], domain_span[1])
    )

    # now extract PDB residues with the correct PDB chain ID, are resolved, have a UniProt crossref and are within the UniProt domain bounds, and do not have "PDB modified", "Conflict" or "Engineered mutation" tags.
    selected_resolved_residues = siftsxml.xpath(
        'entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][not(../residueDetail[contains(text(),"Not_Observed")])]'
        '[../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]][not(../residueDetail[contains(text(),"modified")])]'
        '[not(../residueDetail[contains(text(),"Conflict")])][not(../residueDetail[contains(text(),"mutation")])]' % (chainid, domain_span[0], domain_span[1])
    )

    # all_pdb_domain_residues = siftsxml.xpath(
    #     'entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]]' % (chainid, domain_span[0], domain_span[1])
    # )

    if len(selected_resolved_residues) == 0 or len(selected_residues) == 0:
        return

    # calculate the ratio of resolved residues - if less than a certain amount, discard pdbchain
    ratio_resolved = float(len(selected_resolved_residues)) / float(len(selected_residues))
    if ratio_resolved < ensembler.core.template_acceptable_ratio_resolved_residues:
        return

    # make a single-letter aa code sequence
    full_seq = ''.join([residue.find('../crossRefDb[@dbSource="UniProt"]').get('dbResName') for residue in selected_residues])
    full_pdbresnums = [residue.get('dbResNum') for residue in selected_residues]
    template_seq_resolved = ''.join([residue.find('../crossRefDb[@dbSource="UniProt"]').get('dbResName') for residue in selected_resolved_residues])
    template_pdbresnums_resolved = [residue.get('dbResNum') for residue in selected_resolved_residues]

    # store data
    template_data = TemplateData(
        pdbid=pdbid,
        chainid=chainid,
        templateid=templateid,
        resolved_seq=template_seq_resolved,
        resolved_pdbresnums=template_pdbresnums_resolved,
        full_seq=full_seq,
        full_pdbresnums=full_pdbresnums,
    )

    return template_data


@ensembler.utils.mpirank0only_and_end_with_barrier
def write_template_seqs_to_fasta_file(selected_templates):
    templates_resolved_seqs = [SeqRecord(Seq(template.resolved_seq), id=template.templateid, description=template.templateid) for template in selected_templates]
    templates_full_seqs = [SeqRecord(Seq(template.full_seq), id=template.templateid, description=template.templateid) for template in selected_templates]
    Bio.SeqIO.write(templates_resolved_seqs, os.path.join('templates', 'templates-resolved-seq.fa'), 'fasta')
    Bio.SeqIO.write(templates_full_seqs, os.path.join('templates', 'templates-full-seq.fa'), 'fasta')


@ensembler.utils.mpirank0only_and_end_with_barrier
def extract_template_structures_from_pdb_files(selected_templates, overwrite_structures=False):
    logger.info('Writing template structures...')
    for template in selected_templates:
        pdb_filename = os.path.join(ensembler.core.default_project_dirnames.structures_pdb, template.pdbid + '.pdb.gz')
        template_resolved_filename = os.path.join(ensembler.core.default_project_dirnames.templates_structures_resolved, template.templateid + '.pdb')
        if not overwrite_structures:
            if os.path.exists(template_resolved_filename):
                continue
        ensembler.PDB.extract_residues_by_resnum(template_resolved_filename, pdb_filename, template)


@ensembler.utils.mpirank0only_and_end_with_barrier
def write_gather_templates_from_targetexplorer_metadata(search_string, dbapi_uri, ntemplates, structure_dirs, loopmodel):
    gather_templates_from_targetexplorer_metadata = gen_targetexplorer_metadata(dbapi_uri, search_string)
    gather_templates_from_targetexplorer_metadata['structure_dirs'] = structure_dirs
    gather_templates_metadata = gen_gather_templates_metadata(ntemplates, loopmodel=loopmodel, additional_metadata=gather_templates_from_targetexplorer_metadata)
    project_metadata = ensembler.core.ProjectMetadata(project_stage='gather_templates')
    project_metadata.add_data(gather_templates_metadata)
    project_metadata.write()


@ensembler.utils.mpirank0only_and_end_with_barrier
def write_gather_templates_from_uniprot_metadata(uniprot_query_string, uniprot_domain_regex, ntemplates, structure_dirs, loopmodel):
    gather_templates_from_uniprot_metadata = gen_uniprot_metadata(uniprot_query_string, uniprot_domain_regex)
    gather_templates_from_uniprot_metadata['structure_dirs'] = structure_dirs
    gather_templates_metadata = gen_gather_templates_metadata(ntemplates, loopmodel=loopmodel, additional_metadata=gather_templates_from_uniprot_metadata)
    project_metadata = ensembler.core.ProjectMetadata(project_stage='gather_templates')
    project_metadata.add_data(gather_templates_metadata)
    project_metadata.write()


def gen_gather_templates_metadata(nselected_templates, loopmodel=False, additional_metadata=None):
    if additional_metadata is None:
        additional_metadata = {}
    datestamp = ensembler.core.get_utcnow_formatted()
    metadata = {
        'datestamp': datestamp,
        'ntemplates': str(nselected_templates),
        'loopmodel': loopmodel,
        'python_version': sys.version.split('|')[0].strip(),
        'python_full_version': ensembler.core.literal_str(sys.version),
        'ensembler_version': ensembler.version.short_version,
        'ensembler_commit': ensembler.version.git_revision
    }
    metadata.update(additional_metadata)
    return metadata


def extract_template_pdbchains_from_uniprot_xml(uniprotxml, uniprot_domain_regex, manual_overrides):
    selected_pdbchains = []
    all_uniprot_entries = uniprotxml.findall('entry')
    for entry in all_uniprot_entries:
        entry_name = entry.find('name').text
        if uniprot_domain_regex is not None:
            selected_domains = entry.xpath(
                'feature[@type="domain"][match_regex(@description, "%s")]' % uniprot_domain_regex,
                extensions={(None, 'match_regex'): ensembler.core.xpath_match_regex_case_sensitive}
            )
        else:
            selected_domains = entry.findall('feature[@type="domain"]')

        domain_iter = 0
        for domain in selected_domains:
            domain_id = '%s_D%d' % (entry_name, domain_iter)
            domain_span = [int(domain.find('location/begin').get('position')), int(domain.find('location/end').get('position'))]
            if domain_id in manual_overrides.template.domain_spans:
                domain_span = [int(x) for x in manual_overrides.template.domain_spans[domain_id].split('-')]
            domain_len = domain_span[1] - domain_span[0] + 1
            if manual_overrides.template.min_domain_len is not None and domain_len < manual_overrides.template.min_domain_len:
                continue
            if manual_overrides.template.max_domain_len is not None and domain_len > manual_overrides.template.max_domain_len:
                continue

            domain_iter += 1
            pdbs = domain.getparent().xpath(
                'dbReference[@type="PDB"]/property[@type="method"][@value="X-ray" or @value="NMR"]/..')

            for pdb in pdbs:
                pdbid = pdb.get('id')
                if pdbid in manual_overrides.template.skip_pdbs:
                    continue
                pdb_chain_span_nodes = pdb.findall('property[@type="chains"]')

                for PDB_chain_span_node in pdb_chain_span_nodes:
                    chain_span_string = PDB_chain_span_node.get('value')
                    chain_spans = ensembler.UniProt.parse_uniprot_pdbref_chains(chain_span_string)

                    for chainid in chain_spans.keys():
                        span = chain_spans[chainid]
                        if (span[0] < domain_span[0] + 30) & (span[1] > domain_span[1] - 30):
                            templateid = '%s_%s_%s' % (domain_id, pdbid, chainid)
                            data = {
                                'templateid': templateid,
                                'pdbid': pdbid,
                                'chainid': chainid,
                                'domain_span': domain_span
                            }
                            selected_pdbchains.append(data)
    logger.info('%d PDB chains selected.' % len(selected_pdbchains))
    return selected_pdbchains

structure_type_file_extension_mapper = {'pdb': '.pdb.gz', 'sifts': '.xml.gz'}


def pdbfix_templates(selected_templates, overwrite_structures=False):
    missing_residues_sublist = []
    ntemplates = len(selected_templates)
    for template_index in range(mpistate.rank, ntemplates, mpistate.size):
        missing_residues_sublist.append(pdbfix_template(selected_templates[template_index], overwrite_structures=overwrite_structures))

    missing_residues_gathered = mpistate.comm.gather(missing_residues_sublist, root=0)

    missing_residues_list = []
    if mpistate.rank == 0:
        missing_residues_list = [None] * ntemplates
        for template_index in range(ntemplates):
            missing_residues_list[template_index] = missing_residues_gathered[template_index % mpistate.size][template_index // mpistate.size]

    missing_residues_list = mpistate.comm.bcast(missing_residues_list, root=0)

    return missing_residues_list


def pdbfix_template(template, overwrite_structures=False):
    try:
        log_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-loopmodel-log.yaml')
        if not overwrite_structures:
            if os.path.exists(log_filepath):
                return
        import pdbfixer
        import simtk.openmm.app
        template_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_resolved, template.templateid + '.pdb')
        fixer = pdbfixer.PDBFixer(filename=template_filepath)
        seq_obj = simtk.openmm.app.internal.pdbstructure.Sequence(template.chainid)
        for r in template.full_seq:
            resi3 = Bio.SeqUtils.seq3(r).upper()
            seq_obj.residues.append(resi3)
        fixer.structure.sequences.append(seq_obj)
        fixer.findMissingResidues()
        remove_missing_residues_at_termini(fixer, len_full_seq=len(template.full_seq))
        fixer.findMissingAtoms()
        (newTopology, newPositions, newAtoms, existingAtomMap) = fixer._addAtomsToTopology(True, True)
        fixer.topology = newTopology
        fixer.positions = newPositions
        template_pdbfixed_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-pdbfixed.pdb')
        with open(template_pdbfixed_filepath, 'w') as template_pdbfixed_file:
            simtk.openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, file=template_pdbfixed_file)
        return fixer.missingResidues
    except KeyboardInterrupt:
        raise
    except Exception as e:
        trbk = traceback.format_exc()
        log_filepath = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-pdbfixer-log.yaml'))
        logfile = ensembler.core.LogFile(log_filepath)
        logfile.log({
            'templateid': str(template.templateid),
            'exception': e,   # if Rosetta loopmodel fails, this field should include text sent to stdout and stderr
            'traceback': ensembler.core.literal_str(trbk),
            'mpi_rank': mpistate.rank,
        })
        logger.error('MPI rank %d pdbfixer error for template %s - see logfile' % (mpistate.rank, template.templateid))


def remove_missing_residues_at_termini(fixer, len_full_seq):
    # remove C-terminal missing residues
    if len(fixer.missingResidues) == 0:
        return None
    sorted_missing_residues_keys = sorted(fixer.missingResidues, key=lambda x: x[1])
    last_missing_residues_key = sorted_missing_residues_keys[-1]
    last_missing_residues_start_index = last_missing_residues_key[1]
    last_missing_residues = fixer.missingResidues[last_missing_residues_key]
    nmissing_residues_up_to_last = sum([len(fixer.missingResidues[key]) for key in sorted_missing_residues_keys[:-1]])

    if last_missing_residues_start_index + nmissing_residues_up_to_last + len(last_missing_residues) == len_full_seq:
        fixer.missingResidues.pop(last_missing_residues_key)

    # remove N-terminal missing residues
    fixer.missingResidues.pop((0, 0), None)


def loopmodel_templates(selected_templates, missing_residues, overwrite_structures=False):
    for template_index in range(mpistate.rank, len(selected_templates), mpistate.size):
        template = selected_templates[template_index]
        if mpistate.size > 1:
            logger.info('MPI rank %d modeling missing loops for template %s' % (mpistate.rank, template.templateid))
        else:
            logger.info('Modeling missing loops for template %s' % template.templateid)
        loopmodel_template(template, missing_residues[template_index], overwrite_structures=overwrite_structures)


def loopmodel_template(template, missing_residues, overwrite_structures=False):
    template_filepath = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-pdbfixed.pdb'))
    output_pdb_filepath = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '.pdb'))
    loop_filepath = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '.loop'))
    output_score_filepath = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-loopmodel-score.sc'))
    log_filepath = os.path.abspath(os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-loopmodel-log.yaml'))
    if not overwrite_structures:
        if os.path.exists(log_filepath):
            return
    logfile = ensembler.core.LogFile(log_filepath)
    write_loop_file(template, missing_residues)
    starttime = datetime.datetime.utcnow()
    if len(missing_residues) == 0:
        loopmodel_output = LoopmodelOutput(successful=True, no_missing_residues=True)
    else:
        loopmodel_output = run_loopmodel(template_filepath, loop_filepath, output_pdb_filepath, output_score_filepath)
    if not loopmodel_output.successful:
        logger.error('MPI rank %d Loopmodel error for template %s - see logfile' % (mpistate.rank, template.templateid))
    timedelta = datetime.datetime.utcnow() - starttime
    logfile.log({
        'templateid': str(template.templateid),
        'no_missing_residues': loopmodel_output.no_missing_residues,
        'loopmodel_output': loopmodel_output.output_text,
        'mpi_rank': mpistate.rank,
        'successful': loopmodel_output.successful,
        'exception': loopmodel_output.exception,
        'loopmodel_exception': loopmodel_output.loopmodel_exception,
        'traceback': loopmodel_output.traceback,
        'timing': ensembler.core.strf_timedelta(timedelta),
        })


def write_loop_file(template, missing_residues):
    loop_file_text = ''
    loop_residues_added = 0
    loop_residues_data = [(key[1], len(residues)) for key, residues in missing_residues.iteritems()]
    loop_residues_data = sorted(loop_residues_data, key=lambda x: x[0])
    for loop_residue_data in loop_residues_data:
        residue_number, nresidues = loop_residue_data
        loop_begin = residue_number + loop_residues_added   # 1-based, one residue before the loop
        loop_end = residue_number + nresidues + loop_residues_added + 1   # 1-based, one residue after the loop
        loop_residues_added += nresidues
        # Note that missing residues at termini (which cannot be modeled by Rosetta loopmodel) have already been removed from the PDBFixer.missingResidues dictionary
        loop_file_text += 'LOOP%4d%4d - - 1\n' % (loop_begin, loop_end)
    loop_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '.loop')
    with open(loop_filepath, 'w') as loop_file:
        loop_file.write(loop_file_text)


def run_loopmodel(input_template_pdb_filepath, loop_filepath, output_pdb_filepath, output_score_filepath, loopmodel_executable_filepath=ensembler.core.find_loopmodel_executable(), nmodels_to_build=1):
    temp_dir = tempfile.mkdtemp()
    temp_template_filepath = os.path.join(temp_dir, 'template.pdb')
    temp_loop_filepath = os.path.join(temp_dir, 'template.loop')
    temp_output_model_filepath = os.path.join(temp_dir, 'template_0001.pdb')
    temp_output_score_filepath = os.path.join(temp_dir, 'score.sc')
    minirosetta_database_path = os.environ.get('MINIROSETTA_DATABASE')
    shutil.copy(input_template_pdb_filepath, temp_template_filepath)
    shutil.copy(loop_filepath, temp_loop_filepath)
    try:
        output_text = subprocess.check_output(
            [
                loopmodel_executable_filepath,
                '-database', minirosetta_database_path,
                '-in::file::s', temp_template_filepath,
                '-loops:loop_file', temp_loop_filepath,
                '-out:path:all', temp_dir,
                '-loops:remodel', 'perturb_kic',
                '-loops:refine', 'refine_kic',
                '-ex1',
                '-ex2',
                '-nstruct', '%d' % nmodels_to_build,
                '-loops:max_kic_build_attempts', '100',
                '-in:file:fullatom',
                '-overwrite',
                ],
            stderr=subprocess.STDOUT
        )
        if os.path.exists(temp_output_model_filepath):
            shutil.copy(temp_output_model_filepath, output_pdb_filepath)
            shutil.copy(temp_output_score_filepath, output_score_filepath)
            shutil.rmtree(temp_dir)
            return LoopmodelOutput(output_text=output_text, successful=True)
        else:
            shutil.rmtree(temp_dir)
            return LoopmodelOutput(output_text=output_text, successful=False)
    except KeyboardInterrupt:
        shutil.rmtree(temp_dir)
        raise
    except subprocess.CalledProcessError as e:
        shutil.rmtree(temp_dir)
        return LoopmodelOutput(loopmodel_exception=e.output, trbk=traceback.format_exc(), successful=False)
    except Exception as e:
        shutil.rmtree(temp_dir)
        return LoopmodelOutput(output_text=output_text, exception=e, trbk=traceback.format_exc(), successful=False)


def check_loopmodel_complete_and_successful(template):
    output_pdb_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '.pdb')
    log_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-loopmodel-log.yaml')
    if os.path.exists(log_filepath) and os.path.exists(output_pdb_filepath):
        with open(log_filepath) as log_file:
            log_data = yaml.load(log_file, Loader=ensembler.core.YamlLoader)
            if log_data.get('successful') == True:
                print template.templateid, 'Already pdbfixed'
                return True
    else:
        return False