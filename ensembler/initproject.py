import gzip
import json
import sys
import os
import re

from lxml import etree
import Bio.SeqUtils
import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import ensembler.version

import ensembler
import ensembler.TargetExplorer
import ensembler.uniprot
import ensembler.pdb
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
    def __init__(self, dbapi_uri, search_string='', loglevel=None, run_main=True):
        ensembler.utils.set_loglevel(loglevel)
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
    """Gathers target data from UniProt.

    Attributes
    ----------
    targets: list of Bio.SeqRecord.SeqRecord
        Target ids and sequences
    residue spans: list of [start, end]
        0-based numbering in coordinates of UniProt canonical sequence
    uniprotxml: lxml.etree.Element
        XML returned from UniProt query
    """
    def __init__(self, uniprot_query_string, uniprot_domain_regex=None, save_uniprot_xml=False, loglevel=None, run_main=True):
        ensembler.utils.set_loglevel(loglevel)
        super(GatherTargetsFromUniProt, self).__init__()
        self.uniprot_query_string = uniprot_query_string
        self.uniprot_domain_regex = uniprot_domain_regex
        self._save_uniprot_xml = save_uniprot_xml
        if run_main:
            self._gather_targets()

    @ensembler.utils.notify_when_done
    def _gather_targets(self, write_output_files=True):
        logger.info('Querying UniProt web server...')

        get_uniprot_xml_args = {}
        if self._save_uniprot_xml:
            get_uniprot_xml_args['write_to_filepath'] = 'targets-uniprot.xml'

        self.uniprotxml = ensembler.uniprot.get_uniprot_xml(self.uniprot_query_string, **get_uniprot_xml_args)

        logger.info('Number of entries returned from initial UniProt search: %r\n' % len(self.uniprotxml))
        log_unique_domain_names(self.uniprot_query_string, self.uniprotxml)
        if self.uniprot_domain_regex:
            log_unique_domain_names_selected_by_regex(self.uniprot_domain_regex, self.uniprotxml)
        fasta_ofilepath = os.path.join(ensembler.core.default_project_dirnames.targets, 'targets.fa')
        self._extract_targets_from_uniprot_xml()
        if write_output_files:
            Bio.SeqIO.write(self.targets, fasta_ofilepath, 'fasta')
            self._write_metadata()

    def _extract_targets_from_uniprot_xml(self):
        targets = []
        residue_spans = []
        domain_descriptions = []
        for entry in self.uniprotxml.findall('entry'):
            entry_name = entry.find('name').text
            fullseq = ensembler.core.sequnwrap(entry.find('sequence').text)
            if self.uniprot_domain_regex:
                selected_domains = entry.xpath(
                    'feature[@type="domain"][match_regex(@description, "%s")]' % self.uniprot_domain_regex,
                    extensions={(None, 'match_regex'): ensembler.core.xpath_match_regex_case_sensitive}
                )

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
                    residue_spans.append([start, end])
                    domain_descriptions.append(domain.get('description'))
                    domain_iter += 1

            else:
                targetid = entry_name
                targets.append(SeqRecord(Seq(fullseq), id=targetid, description=targetid))
                residue_spans.append([0, len(fullseq)-1])

        self.targets = targets
        self.residue_spans = residue_spans
        self.domain_descriptions = domain_descriptions

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


def gen_pdb_metadata(pdbids, uniprot_domain_regex, chainids):
    metadata = {
        'method': 'PDB',
        'gather_from_pdb': {
            'pdbids': pdbids,
            'uniprot_domain_regex': uniprot_domain_regex,
            'chainids': chainids,
        }
    }
    return metadata


def log_unique_domain_names(uniprot_query_string, uniprotxml):
    # Example query string: 'domain:"Protein kinase" AND reviewed:yes'
    domain_match = re.search('domain:([\"\'].*[\"\'])', uniprot_query_string)
    if domain_match and len(domain_match.groups()) > 0:
        query_string_domain_selection = domain_match.groups()[0].replace('\'', '').replace('\"', '')
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
def gather_templates_from_targetexplorer(dbapi_uri, search_string='', structure_dirs=None, loglevel=None):
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
    ensembler.utils.set_loglevel(loglevel)
    manual_overrides = ensembler.core.ManualOverrides()
    templates_json = get_targetexplorer_templates_json(dbapi_uri, search_string)
    selected_pdbchains = extract_template_pdbchains_from_targetexplorer_json(templates_json, manual_overrides=manual_overrides)
    for pdbchain in selected_pdbchains:
        get_pdb_and_sifts_files(pdbchain['pdbid'], structure_dirs)

    selected_templates = extract_template_pdb_chain_residues(selected_pdbchains)
    write_template_seqs_to_fasta_file(selected_templates)
    extract_template_structures_from_pdb_files(selected_templates)
    write_gather_templates_from_targetexplorer_metadata(search_string, dbapi_uri, len(selected_templates), structure_dirs)


@ensembler.utils.notify_when_done
def gather_templates_from_uniprot(uniprot_query_string, uniprot_domain_regex=None, structure_dirs=None, pdbids=None, chainids=None, loglevel=None):
    """# Searches UniProt for a set of template proteins with a user-defined
    query string, then saves IDs, sequences and structures."""
    ensembler.utils.set_loglevel(loglevel)
    manual_overrides = ensembler.core.ManualOverrides()
    selected_pdbchains = None
    if mpistate.rank == 0:
        uniprotxml = ensembler.uniprot.get_uniprot_xml(uniprot_query_string)
        log_unique_domain_names(uniprot_query_string, uniprotxml)
        if uniprot_domain_regex is not None:
            log_unique_domain_names_selected_by_regex(uniprot_domain_regex, uniprotxml)

        selected_pdbchains = extract_template_pdbchains_from_uniprot_xml(uniprotxml, uniprot_domain_regex=uniprot_domain_regex, manual_overrides=manual_overrides, specified_pdbids=pdbids, specified_chainids=chainids)
        for pdbchain in selected_pdbchains:
            get_pdb_and_sifts_files(pdbchain['pdbid'], structure_dirs)

    selected_pdbchains = mpistate.comm.bcast(selected_pdbchains, root=0)
    logger.debug('Selected PDB chains: {0}'.format([pdbchain['templateid'] for pdbchain in selected_pdbchains]))

    selected_templates = extract_template_pdb_chain_residues(selected_pdbchains)
    write_template_seqs_to_fasta_file(selected_templates)
    extract_template_structures_from_pdb_files(selected_templates)
    write_gather_templates_from_uniprot_metadata(uniprot_query_string, uniprot_domain_regex, len(selected_templates), structure_dirs)


@ensembler.utils.notify_when_done
def gather_templates_from_pdb(pdbids, uniprot_domain_regex=None, chainids=None, structure_dirs=None, loglevel=None):
    """
    :param pdbids: list of str
    :param uniprot_domain_regex: str
    :param chainids: dict {pdbid (str): [chainid (str)]}
    :param structure_dirs: list of str
    :return:
    """
    ensembler.utils.set_loglevel(loglevel)
    manual_overrides = ensembler.core.ManualOverrides()
    selected_pdbchains = None
    if mpistate.rank == 0:
        for pdbid in pdbids:
            get_pdb_and_sifts_files(pdbid, structure_dirs)
        uniprot_acs = extract_uniprot_acs_from_sifts_files(pdbids)
        logger.debug('Extracted UniProt ACs: {0}'.format(uniprot_acs))
        uniprot_ac_query_string = ensembler.uniprot.build_uniprot_query_string_from_acs(uniprot_acs)
        uniprotxml = ensembler.uniprot.get_uniprot_xml(uniprot_ac_query_string)
        selected_pdbchains = extract_template_pdbchains_from_uniprot_xml(uniprotxml, uniprot_domain_regex=uniprot_domain_regex, manual_overrides=manual_overrides, specified_pdbids=pdbids, specified_chainids=chainids)

    selected_pdbchains = mpistate.comm.bcast(selected_pdbchains, root=0)
    logger.debug('Selected PDB chains: {0}'.format([pdbchain['templateid'] for pdbchain in selected_pdbchains]))

    selected_templates = extract_template_pdb_chain_residues(selected_pdbchains)
    write_template_seqs_to_fasta_file(selected_templates)
    extract_template_structures_from_pdb_files(selected_templates)
    write_gather_templates_from_pdb_metadata(pdbids, uniprot_domain_regex, len(selected_templates), chainids, structure_dirs)


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
    pdbgz_page = ensembler.pdb.retrieve_pdb(pdbid, compressed='yes')
    with open(project_pdb_filepath, 'w') as pdbgz_file:
        pdbgz_file.write(pdbgz_page)


def download_sifts_file(pdbid, project_sifts_filepath):
    logger.info('Downloading sifts file for: %s', pdbid)
    sifts_page = ensembler.pdb.retrieve_sifts(pdbid)
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
    residue_span = pdbchain['residue_span']   # UniProt coords

    sifts_filepath = os.path.join('structures', 'sifts', pdbid + '.xml.gz')
    siftsxml = parse_sifts_xml(sifts_filepath)

    add_pdb_modified_xml_tags_to_residues(siftsxml)

    # domain_span_sifts_coords = [
    #     int(siftsxml.find('entity/segment/listResidue/residue/crossRefDb[@dbSource="UniProt"][@dbChainId="%s"][@dbResNum="%d"]/..' % (chainid, residue_span[0])).get('dbResNum')),
    #     int(siftsxml.find('entity/segment/listResidue/residue/crossRefDb[@dbSource="UniProt"][@dbChainId="%s"][@dbResNum="%d"]/..' % (chainid, residue_span[1])).get('dbResNum')),
    # ]

    # TODO once working, check whether this fills in all loops
    # An alternative approach would be to just take the UniProt sequence specified by the domain span
    selected_residues = siftsxml.xpath(
        'entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"]'
        '[../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]]' % (chainid, residue_span[0], residue_span[1])
    )

    # now extract PDB residues which have the correct PDB chain ID, are resolved, have a UniProt crossref and are within the UniProt domain bounds, and do not have "PDB modified", "Conflict" or "Engineered mutation" tags.
    selected_resolved_residues = siftsxml.xpath(
        'entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][not(../residueDetail[contains(text(),"Not_Observed")])]'
        '[../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]][not(../residueDetail[contains(text(),"modified")])]'
        '[not(../residueDetail[contains(text(),"Conflict")])][not(../residueDetail[contains(text(),"mutation")])]' % (chainid, residue_span[0], residue_span[1])
    )

    # second stage of filtering to remove residues which conflict with the UniProt resname, but are not annotated as such
    selected_resolved_residues = [r for r in selected_resolved_residues if Bio.SeqUtils.seq1(r.get('dbResName')) == r.find('../crossRefDb[@dbSource="UniProt"]').get('dbResName')]

    # all_pdb_domain_residues = siftsxml.xpath(
    #     'entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]]' % (chainid, residue_span[0], residue_span[1])
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
def extract_template_structures_from_pdb_files(selected_templates):
    logger.info('Writing template structures...')
    for template in selected_templates:
        pdb_filename = os.path.join(ensembler.core.default_project_dirnames.structures_pdb, template.pdbid + '.pdb.gz')
        template_resolved_filename = os.path.join(ensembler.core.default_project_dirnames.templates_structures_resolved, template.templateid + '.pdb')
        ensembler.pdb.extract_residues_by_resnum(template_resolved_filename, pdb_filename, template)


@ensembler.utils.mpirank0only_and_end_with_barrier
def write_gather_templates_from_targetexplorer_metadata(search_string, dbapi_uri, ntemplates, structure_dirs):
    gather_templates_from_targetexplorer_metadata = gen_targetexplorer_metadata(dbapi_uri, search_string)
    gather_templates_from_targetexplorer_metadata['structure_dirs'] = structure_dirs
    gather_templates_metadata = gen_gather_templates_metadata(ntemplates, additional_metadata=gather_templates_from_targetexplorer_metadata)
    project_metadata = ensembler.core.ProjectMetadata(project_stage='gather_templates')
    project_metadata.add_data(gather_templates_metadata)
    project_metadata.write()


@ensembler.utils.mpirank0only_and_end_with_barrier
def write_gather_templates_from_uniprot_metadata(uniprot_query_string, uniprot_domain_regex, ntemplates, structure_dirs):
    gather_templates_from_uniprot_metadata = gen_uniprot_metadata(uniprot_query_string, uniprot_domain_regex)
    gather_templates_from_uniprot_metadata['structure_dirs'] = structure_dirs
    gather_templates_metadata = gen_gather_templates_metadata(ntemplates, additional_metadata=gather_templates_from_uniprot_metadata)
    project_metadata = ensembler.core.ProjectMetadata(project_stage='gather_templates')
    project_metadata.add_data(gather_templates_metadata)
    project_metadata.write()


@ensembler.utils.mpirank0only_and_end_with_barrier
def write_gather_templates_from_pdb_metadata(pdbids, uniprot_domain_regex, ntemplates, chainids, structure_dirs):
    gather_templates_from_pdb_metadata = gen_pdb_metadata(pdbids, uniprot_domain_regex, chainids)
    gather_templates_from_pdb_metadata['structure_dirs'] = structure_dirs
    gather_templates_metadata = gen_gather_templates_metadata(ntemplates, additional_metadata=gather_templates_from_pdb_metadata)
    project_metadata = ensembler.core.ProjectMetadata(project_stage='gather_templates')
    project_metadata.add_data(gather_templates_metadata)
    project_metadata.write()


def gen_gather_templates_metadata(nselected_templates, additional_metadata=None):
    if additional_metadata is None:
        additional_metadata = {}
    datestamp = ensembler.core.get_utcnow_formatted()
    metadata = {
        'datestamp': datestamp,
        'ntemplates': str(nselected_templates),
        'python_version': sys.version.split('|')[0].strip(),
        'python_full_version': ensembler.core.literal_str(sys.version),
        'ensembler_version': ensembler.version.short_version,
        'ensembler_commit': ensembler.version.git_revision
    }
    metadata.update(additional_metadata)
    return metadata


def dep_extract_template_pdbchains_from_uniprot_xml(uniprotxml, uniprot_domain_regex=None, manual_overrides=None, specified_pdbids=None, specified_chainids=None):
    selected_pdbchains = []
    all_uniprot_entries = uniprotxml.findall('entry')
    for entry in all_uniprot_entries:
        entry_name = entry.find('name').text
        if uniprot_domain_regex:
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
            if manual_overrides and domain_id in manual_overrides.template.domain_spans:
                domain_span = [int(x) for x in manual_overrides.template.domain_spans[domain_id].split('-')]
            domain_len = domain_span[1] - domain_span[0] + 1
            if manual_overrides and manual_overrides.template.min_domain_len is not None and domain_len < manual_overrides.template.min_domain_len:
                continue
            if manual_overrides and manual_overrides.template.max_domain_len is not None and domain_len > manual_overrides.template.max_domain_len:
                continue

            domain_iter += 1
            pdbs = domain.getparent().xpath(
                'dbReference[@type="PDB"]/property[@type="method"][@value="X-ray" or @value="NMR"]/..')

            for pdb in pdbs:
                pdbid = pdb.get('id')
                if manual_overrides and pdbid in manual_overrides.template.skip_pdbs:
                    continue
                if specified_pdbids and pdbid not in specified_pdbids:
                    continue
                pdb_chain_span_nodes = pdb.findall('property[@type="chains"]')

                for PDB_chain_span_node in pdb_chain_span_nodes:
                    chain_span_string = PDB_chain_span_node.get('value')
                    chain_spans = ensembler.uniprot.parse_uniprot_pdbref_chains(chain_span_string)

                    for chainid in chain_spans.keys():
                        if specified_chainids and len(specified_chainids[pdbid]) > 0 and chainid not in specified_chainids[pdbid]:
                            continue
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


def extract_template_pdbchains_from_uniprot_xml(uniprotxml, uniprot_domain_regex=None, manual_overrides=None, specified_pdbids=None, specified_chainids=None):
    selected_pdbchains = []
    all_uniprot_entries = uniprotxml.findall('entry')
    for entry in all_uniprot_entries:
        entry_name = entry.find('name').text
        if uniprot_domain_regex:
            selected_domains = entry.xpath(
                'feature[@type="domain"][match_regex(@description, "%s")]' % uniprot_domain_regex,
                extensions={(None, 'match_regex'): ensembler.core.xpath_match_regex_case_sensitive}
            )

            domain_iter = 0
            for domain in selected_domains:
                domain_id = '%s_D%d' % (entry_name, domain_iter)
                domain_span = [int(domain.find('location/begin').get('position')), int(domain.find('location/end').get('position'))]
                if manual_overrides and domain_id in manual_overrides.template.domain_spans:
                    domain_span = [int(x) for x in manual_overrides.template.domain_spans[domain_id].split('-')]
                domain_len = domain_span[1] - domain_span[0] + 1
                if manual_overrides and manual_overrides.template.min_domain_len is not None and domain_len < manual_overrides.template.min_domain_len:
                    continue
                if manual_overrides and manual_overrides.template.max_domain_len is not None and domain_len > manual_overrides.template.max_domain_len:
                    continue

                domain_iter += 1
                pdbs = domain.getparent().xpath(
                    'dbReference[@type="PDB"]/property[@type="method"][@value="X-ray" or @value="NMR"]/..'
                )

                for pdb in pdbs:
                    pdbid = pdb.get('id')
                    if manual_overrides and pdbid in manual_overrides.template.skip_pdbs:
                        continue
                    if specified_pdbids and pdbid not in specified_pdbids:
                        continue
                    pdb_chain_span_nodes = pdb.findall('property[@type="chains"]')

                    for pdb_chain_span_node in pdb_chain_span_nodes:
                        chain_span_string = pdb_chain_span_node.get('value')
                        chain_spans = ensembler.uniprot.parse_uniprot_pdbref_chains(chain_span_string)

                        for chainid in chain_spans.keys():
                            if specified_chainids and len(specified_chainids[pdbid]) > 0 and chainid not in specified_chainids[pdbid]:
                                continue
                            span = chain_spans[chainid]
                            if (span[0] < domain_span[0] + 30) & (span[1] > domain_span[1] - 30):
                                templateid = '%s_%s_%s' % (domain_id, pdbid, chainid)
                                data = {
                                    'templateid': templateid,
                                    'pdbid': pdbid,
                                    'chainid': chainid,
                                    'residue_span': domain_span
                                }
                                selected_pdbchains.append(data)

        else:
            pdbs = entry.xpath(
                'dbReference[@type="PDB"]/property[@type="method"][@value="X-ray" or @value="NMR"]/..'
            )

            for pdb in pdbs:
                pdbid = pdb.get('id')
                if manual_overrides and pdbid in manual_overrides.template.skip_pdbs:
                    continue
                if specified_pdbids and pdbid not in specified_pdbids:
                    continue
                pdb_chain_span_nodes = pdb.findall('property[@type="chains"]')

                for pdb_chain_span_node in pdb_chain_span_nodes:
                    chain_span_string = pdb_chain_span_node.get('value')
                    chain_spans = ensembler.uniprot.parse_uniprot_pdbref_chains(chain_span_string)

                    for chainid in chain_spans.keys():
                        if specified_chainids and len(specified_chainids[pdbid]) > 0 and chainid not in specified_chainids[pdbid]:
                            continue
                        span = chain_spans[chainid]
                        templateid = '%s_%s_%s' % (entry_name, pdbid, chainid)
                        data = {
                            'templateid': templateid,
                            'pdbid': pdbid,
                            'chainid': chainid,
                            'residue_span': span
                        }
                        selected_pdbchains.append(data)

    logger.info('%d PDB chains selected.' % len(selected_pdbchains))
    return selected_pdbchains


structure_type_file_extension_mapper = {'pdb': '.pdb.gz', 'sifts': '.xml.gz'}


def extract_uniprot_acs_from_sifts_files(pdbids):
    uniprot_acs = []
    for pdbid in pdbids:
        sifts_filepath = os.path.join(ensembler.core.default_project_dirnames.structures_sifts, pdbid + '.xml.gz')
        parser = etree.XMLParser(huge_tree=True)
        siftsxml = etree.parse(sifts_filepath, parser).getroot()
        uniprot_acs += ensembler.pdb.extract_uniprot_acs_from_sifts_xml(siftsxml)
    return list(set(uniprot_acs))