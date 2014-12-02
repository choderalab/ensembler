import gzip
import json
import sys
import os
import shutil
import logging
import subprocess
import traceback
import tempfile
from lxml import etree
import msmseeder
import msmseeder.TargetExplorer
import msmseeder.UniProt
import msmseeder.PDB
import msmseeder.version
from msmseeder.core import construct_fasta_str
from msmseeder.utils import file_exists_and_not_empty
from collections import namedtuple
from msmseeder.core import mpistate

logger = logging.getLogger('info')


class TemplateData:
    def __init__(self, pdbid=None, chainid=None, templateid=None, resolved_seq=None, resolved_pdbresnums=None, full_seq=None, full_pdbresnums=None):
        self.pdbid = pdbid
        self.chainid = chainid
        self.templateid = templateid
        self.resolved_seq = resolved_seq
        self.resolved_pdbresnums = resolved_pdbresnums
        self.full_seq = full_seq
        self.full_pdbresnums = full_pdbresnums


@msmseeder.utils.notify_when_done
def initproject(project_toplevel_dir):
    """Initialize MSMSeeder project within the given directory. Creates
    necessary subdirectories and a project metadata .yaml file.
    :param project_toplevel_dir: str
    """
    create_project_dirs(project_toplevel_dir)
    write_init_metadata(project_toplevel_dir)


@msmseeder.utils.notify_when_done
def gather_targets_from_targetexplorer(dbapi_uri, search_string=''):
    """Gather protein target data from a TargetExplorer DB network API.
    Pass the URI for the database API and a search string.
    The search string uses SQLAlchemy syntax and standard TargetExplorer
    frontend data fields.
    Example:
    dbapi_uri='http://plfah2.mskcc.org/kinomeDBAPI'
    search_string='species="Human"'

    To select all domains within the database:
    search_string=''
    """
    manual_overrides = msmseeder.core.ManualOverrides()
    domain_span_overrides_present = True if len(manual_overrides.target.domain_spans) > 0 else False

    targets_json = get_targetexplorer_targets_json(dbapi_uri, search_string, domain_span_overrides_present)
    targets = extract_targets_from_targetexplorer_json(targets_json, manual_overrides=manual_overrides)
    write_seqs_to_fasta_file(targets)

    write_gather_targets_from_targetexplorer_metadata(dbapi_uri, search_string, len(targets))


@msmseeder.utils.notify_when_done
def gather_targets_from_uniprot(uniprot_query_string, uniprot_domain_regex):
    """Searches UniProt for a set of target proteins with a user-defined
    query string, then saves target IDs and sequences."""
    manual_overrides = msmseeder.core.ManualOverrides()
    uniprotxml = get_uniprot_xml(uniprot_query_string)
    log_unique_domain_names(uniprot_query_string, uniprotxml)
    if uniprot_domain_regex != None:
        log_unique_domain_names_selected_by_regex(uniprot_domain_regex, uniprotxml)
    targets = extract_targets_from_uniprot_xml(uniprotxml, uniprot_domain_regex, manual_overrides)
    write_seqs_to_fasta_file(targets)
    write_gather_targets_from_uniprot_metadata(uniprot_query_string, uniprot_domain_regex, len(targets))


def create_project_dirs(project_toplevel_dir):
    msmseeder.utils.create_dir(os.path.join(project_toplevel_dir, msmseeder.core.default_project_dirnames.targets))
    msmseeder.utils.create_dir(os.path.join(project_toplevel_dir, msmseeder.core.default_project_dirnames.templates))
    msmseeder.utils.create_dir(os.path.join(project_toplevel_dir, msmseeder.core.default_project_dirnames.structures))
    msmseeder.utils.create_dir(os.path.join(project_toplevel_dir, msmseeder.core.default_project_dirnames.models))
    msmseeder.utils.create_dir(os.path.join(project_toplevel_dir, msmseeder.core.default_project_dirnames.packaged_models))
    msmseeder.utils.create_dir(os.path.join(project_toplevel_dir, msmseeder.core.default_project_dirnames.structures_pdb))
    msmseeder.utils.create_dir(os.path.join(project_toplevel_dir, msmseeder.core.default_project_dirnames.structures_sifts))
    msmseeder.utils.create_dir(os.path.join(project_toplevel_dir, msmseeder.core.default_project_dirnames.templates_structures_resolved))
    msmseeder.utils.create_dir(os.path.join(project_toplevel_dir, msmseeder.core.default_project_dirnames.templates_structures_modeled_loops))


def write_init_metadata(project_toplevel_dir):
    project_metadata = msmseeder.core.ProjectMetadata(project_stage='init', project_toplevel_dir=project_toplevel_dir)
    init_metadata = gen_init_metadata(project_toplevel_dir)
    project_metadata.add_data(init_metadata)
    project_metadata.write()


def gen_init_metadata(project_toplevel_dir):
    datestamp = msmseeder.core.get_utcnow_formatted()
    metadata_dict = {
        'datestamp': datestamp,
        'init_path': os.path.abspath(project_toplevel_dir),
        'python_version': sys.version.split('|')[0].strip(),
        'python_full_version': msmseeder.core.literal_str(sys.version),
        'msmseeder_version': msmseeder.version.short_version,
        'msmseeder_commit': msmseeder.version.git_revision
    }
    return metadata_dict


def get_targetexplorer_targets_json(dbapi_uri, search_string, domain_span_overrides_present=False):
    """
    :param dbapi_uri: str
    :param search_string: str
    :param domain_span_overrides_present: bool
    :return: list containing nested lists and dicts
    """
    if domain_span_overrides_present:
        targetexplorer_jsonstr = msmseeder.TargetExplorer.query_targetexplorer(
            dbapi_uri, search_string, return_data='domain_seqs,seqs'
        )
    else:
        targetexplorer_jsonstr = msmseeder.TargetExplorer.query_targetexplorer(
            dbapi_uri, search_string, return_data='domain_seqs'
        )
    targetexplorer_json = json.loads(targetexplorer_jsonstr)
    return targetexplorer_json


def get_uniprot_xml(uniprot_query_string):
    logger.info('Querying UniProt web server...')
    uniprotxmlstring = msmseeder.UniProt.retrieve_uniprot(uniprot_query_string)
    parser = etree.XMLParser(huge_tree=True)
    uniprotxml = etree.fromstring(uniprotxmlstring, parser)
    logger.info('Number of entries returned from initial UniProt search: %r\n' % len(uniprotxml))
    return uniprotxml


def extract_targets_from_targetexplorer_json(targets_json, manual_overrides=msmseeder.core.ManualOverrides()):
    """
    :param targets_json:
    :param manual_overrides:
    :return: list of tuples: [(id, seq), (id, seq), ...]
    """
    targets = []
    for target in targets_json['results']:
        for target_domain in target['domains']:
            targetid = target_domain.get('targetid')
            targetseq = target_domain.get('sequence')
            # domain span override
            if targetid in manual_overrides.target.domain_spans:
                fullseq = target.get('sequence')
                start, end = [int(x) - 1 for x in manual_overrides.target.domain_spans[targetid].split('-')]
                targetseq = fullseq[start:end + 1]
            targets.append((targetid, targetseq))

    return targets


def extract_targets_from_uniprot_xml(uniprotxml, uniprot_domain_regex, manual_overrides):
    targets = []
    for entry in uniprotxml.findall('entry'):
        entry_name = entry.find('name').text
        fullseq = msmseeder.core.sequnwrap(entry.find('sequence').text)
        if uniprot_domain_regex is not None:
            selected_domains = entry.xpath(
                'feature[@type="domain"][match_regex(@description, "%s")]' % uniprot_domain_regex,
                extensions={(None, 'match_regex'): msmseeder.core.xpath_match_regex_case_sensitive}
            )
        else:
            selected_domains = entry.findall('feature[@type="domain"]')

        domain_iter = 0
        for domain in selected_domains:
            targetid = '%s_D%d' % (entry_name, domain_iter)
            # domain span override
            if targetid in manual_overrides.target.domain_spans:
                start, end = [int(x) - 1 for x in manual_overrides.target.domain_spans[targetid].split('-')]
            else:
                start, end = [int(domain.find('location/begin').get('position')) - 1,
                              int(domain.find('location/end').get('position')) - 1]
            targetseq = fullseq[start:end + 1]
            targets.append((targetid, targetseq))
            domain_iter += 1

    return targets


def write_seqs_to_fasta_file(targets, fasta_ofilepath=os.path.join('targets', 'targets.fa')):
    """
    :param targets: list of tuples [(id, seq), (id, seq), ...]
    :param fasta_ofilepath: str
    """
    logger.info('Writing target data to FASTA file "%s"...' % fasta_ofilepath)
    with open(fasta_ofilepath, 'w') as fasta_ofile:
        for target in targets:
            targetid, targetseq = target
            target_fasta_string = construct_fasta_str(targetid, targetseq)
            fasta_ofile.write(target_fasta_string)


def gen_metadata_gather_from_targetexplorer(search_string, dbapi_uri):
    db_metadata = get_targetexplorer_db_metadata(dbapi_uri)
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


def gen_metadata_gather_from_uniprot(uniprot_query_string, uniprot_domain_regex):
    metadata = {
        'method': 'UniProt',
        'gather_from_uniprot': {
            'uniprot_query_string': uniprot_query_string,
            'uniprot_domain_regex': uniprot_domain_regex,
        }
    }
    return metadata


def get_targetexplorer_db_metadata(dbapi_uri):
    db_metadata_jsonstr = msmseeder.TargetExplorer.get_targetexplorer_metadata(dbapi_uri)
    return json.loads(db_metadata_jsonstr)


def write_gather_targets_from_targetexplorer_metadata(dbapi_uri, search_string, ntargets):
    gather_from_targetexplorer_metadata = gen_metadata_gather_from_targetexplorer(search_string, dbapi_uri)
    gather_targets_metadata = gen_gather_targets_metadata(ntargets, additional_metadata=gather_from_targetexplorer_metadata)
    project_metadata = msmseeder.core.ProjectMetadata(project_stage='gather_targets')
    project_metadata.add_data(gather_targets_metadata)
    project_metadata.write()


def write_gather_targets_from_uniprot_metadata(uniprot_query_string, uniprot_domain_regex, ntargets):
    gather_from_uniprot_metadata = gen_metadata_gather_from_uniprot(uniprot_query_string, uniprot_domain_regex)
    gather_targets_metadata = gen_gather_targets_metadata(ntargets, additional_metadata=gather_from_uniprot_metadata)
    project_metadata = msmseeder.core.ProjectMetadata(project_stage='gather_targets')
    project_metadata.add_data(gather_targets_metadata)
    project_metadata.write()


def gen_gather_targets_metadata(ntarget_domains, additional_metadata=None):
    if additional_metadata is None:
        additional_metadata = {}
    datestamp = msmseeder.core.get_utcnow_formatted()
    metadata = {
        'datestamp': datestamp,
        'ntargets': str(ntarget_domains),
        'python_version': sys.version.split('|')[0].strip(),
        'python_full_version': msmseeder.core.literal_str(sys.version),
        'msmseeder_version': msmseeder.version.short_version,
        'msmseeder_commit': msmseeder.version.git_revision,
    }
    metadata.update(additional_metadata)
    return metadata

@msmseeder.utils.mpirank0only_and_end_with_barrier
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
                (None, 'match_regex'): msmseeder.core.xpath_match_regex_case_insensitive
            }
        )
        uniprot_unique_domain_names = set([domain.get('description') for domain in uniprot_query_string_domains])
        logger.info('Set of unique domain names selected by the domain selector \'%s\'during the initial UniProt search:\n%s\n'
              % (query_string_domain_selection, uniprot_unique_domain_names))

    else:
        uniprot_domains = uniprotxml.xpath('entry/feature[@type="domain"]')
        uniprot_unique_domain_names = set([domain.get('description') for domain in uniprot_domains])
        logger.info('Set of unique domain names returned from the initial UniProt search using the query string \'%s\':\n%s\n'
              % (uniprot_query_string, uniprot_unique_domain_names))


@msmseeder.utils.mpirank0only_and_end_with_barrier
def log_unique_domain_names_selected_by_regex(uniprot_domain_regex, uniprotxml):
    regex_matched_domains = uniprotxml.xpath(
        'entry/feature[@type="domain"][match_regex(@description, "%s")]' % uniprot_domain_regex,
        extensions={(None, 'match_regex'): msmseeder.core.xpath_match_regex_case_sensitive}
    )
    regex_matched_domains_unique_names = set([domain.get('description') for domain in regex_matched_domains])
    logger.info('Unique domain names selected after searching with the case-sensitive regex string \'%s\':\n%s\n'
        % (uniprot_domain_regex, regex_matched_domains_unique_names))


@msmseeder.utils.notify_when_done
def gather_templates_from_targetexplorer(dbapi_uri, search_string='', structure_dirs=None, loopmodel=True):
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
    manual_overrides = msmseeder.core.ManualOverrides()
    templates_json = get_targetexplorer_templates_json(dbapi_uri, search_string)
    selected_pdbchains = extract_template_pdbchains_from_targetexplorer_json(templates_json, manual_overrides=manual_overrides)
    for pdbchain in selected_pdbchains:
        get_pdb_and_sifts_files(pdbchain['pdbid'], structure_dirs)

    selected_templates = extract_template_pdb_chain_residues(selected_pdbchains)
    write_template_seqs_to_fasta_file(selected_templates)
    extract_template_structures_from_pdb_files(selected_templates)
    if loopmodel:
        missing_residues = pdbfix_templates(selected_templates)
        loopmodel_templates(selected_templates, missing_residues)
    write_gather_templates_from_targetexplorer_metadata(search_string, dbapi_uri, len(selected_templates), structure_dirs, loopmodel)


@msmseeder.utils.notify_when_done
def gather_templates_from_uniprot(uniprot_query_string, uniprot_domain_regex, structure_dirs=None, loopmodel=True):
    """# Searches UniProt for a set of template proteins with a user-defined
    query string, then saves IDs, sequences and structures."""
    manual_overrides = msmseeder.core.ManualOverrides()
    selected_pdbchains = None
    if mpistate.rank == 0:
        uniprotxml = get_uniprot_xml(uniprot_query_string)
        log_unique_domain_names(uniprot_query_string, uniprotxml)
        if uniprot_domain_regex is not None:
            log_unique_domain_names_selected_by_regex(uniprot_domain_regex, uniprotxml)

        selected_pdbchains = extract_template_pdbchains_from_uniprot_xml(uniprotxml, uniprot_domain_regex, manual_overrides)
        for pdbchain in selected_pdbchains:
            get_pdb_and_sifts_files(pdbchain['pdbid'], structure_dirs)

    selected_pdbchains = mpistate.comm.bcast(selected_pdbchains, root=0)

    selected_templates = extract_template_pdb_chain_residues(selected_pdbchains)
    write_template_seqs_to_fasta_file(selected_templates)
    extract_template_structures_from_pdb_files(selected_templates)
    if loopmodel:
        missing_residues = pdbfix_templates(selected_templates)
        loopmodel_templates(selected_templates, missing_residues)
    write_gather_templates_from_uniprot_metadata(uniprot_query_string, uniprot_domain_regex, len(selected_templates), structure_dirs, loopmodel)


def get_targetexplorer_templates_json(dbapi_uri, search_string):
    """
    :param dbapi_uri: str
    :param search_string: str
    :return: list containing nested lists and dicts
    """
    targetexplorer_json = None
    if mpistate.rank == 0:
        targetexplorer_jsonstr = msmseeder.TargetExplorer.query_targetexplorer(
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
    pdbgz_page = msmseeder.PDB.retrieve_pdb(pdbid, compressed='yes')
    with open(project_pdb_filepath, 'w') as pdbgz_file:
        pdbgz_file.write(pdbgz_page)


def download_sifts_file(pdbid, project_sifts_filepath):
    logger.info('Downloading sifts file for: %s', pdbid)
    sifts_page = msmseeder.PDB.retrieve_sifts(pdbid)
    with gzip.open(project_sifts_filepath, 'wb') as project_sifts_file:
        project_sifts_file.write(sifts_page)


@msmseeder.utils.mpirank0only_and_end_with_barrier
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

    # now extract PDB residues with the correct PDB chain ID, are observed, have a UniProt crossref and are within the UniProt domain bounds, and do not have "PDB modified", "Conflict" or "Engineered mutation" tags.
    selected_observed_residues = siftsxml.xpath(
        'entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][not(../residueDetail[contains(text(),"Not_Observed")])]'
        '[../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]][not(../residueDetail[contains(text(),"modified")])]'
        '[not(../residueDetail[contains(text(),"Conflict")])][not(../residueDetail[contains(text(),"mutation")])]' % (chainid, domain_span[0], domain_span[1])
    )

    # all_pdb_domain_residues = siftsxml.xpath(
    #     'entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]]' % (chainid, domain_span[0], domain_span[1])
    # )

    if len(selected_observed_residues) == 0 or len(selected_residues) == 0:
        return

    # calculate the ratio of observed residues - if less than a certain amount, discard pdbchain
    ratio_observed = float(len(selected_observed_residues)) / float(len(selected_residues))
    if ratio_observed < msmseeder.core.template_acceptable_ratio_observed_residues:
        return

    # make a single-letter aa code sequence
    full_seq = ''.join([residue.find('../crossRefDb[@dbSource="UniProt"]').get('dbResName') for residue in selected_residues])
    full_pdbresnums = [residue.get('dbResNum') for residue in selected_residues]
    template_seq_observed = ''.join([residue.find('../crossRefDb[@dbSource="UniProt"]').get('dbResName') for residue in selected_observed_residues])
    template_pdbresnums_observed = [residue.get('dbResNum') for residue in selected_observed_residues]

    # store data
    template_data = TemplateData(
        pdbid=pdbid,
        chainid=chainid,
        templateid=templateid,
        resolved_seq=template_seq_observed,
        resolved_pdbresnums=template_pdbresnums_observed,
        full_seq=full_seq,
        full_pdbresnums=full_pdbresnums,
    )

    return template_data


@msmseeder.utils.mpirank0only_and_end_with_barrier
def write_template_seqs_to_fasta_file(selected_templates):
    selected_template_seq_tuples = [(template.templateid, template.full_seq) for template in selected_templates]
    selected_template_seq_observed_tuples = [(template.templateid, template.resolved_seq) for template in selected_templates]
    write_seqs_to_fasta_file(selected_template_seq_tuples, fasta_ofilepath=os.path.join('templates', 'templates-full-seq.fa'))
    write_seqs_to_fasta_file(selected_template_seq_observed_tuples, fasta_ofilepath=os.path.join('templates', 'templates-resolved.fa'))


@msmseeder.utils.mpirank0only_and_end_with_barrier
def extract_template_structures_from_pdb_files(selected_templates):
    logger.info('Writing template structures...')
    for template in selected_templates:
        pdb_filename = os.path.join(msmseeder.core.default_project_dirnames.structures_pdb, template.pdbid + '.pdb.gz')
        template_observed_filename = os.path.join(msmseeder.core.default_project_dirnames.templates_structures_resolved, template.templateid + '.pdb')
        msmseeder.PDB.extract_residues_by_resnum(template_observed_filename, pdb_filename, template)


@msmseeder.utils.mpirank0only_and_end_with_barrier
def write_gather_templates_from_targetexplorer_metadata(search_string, dbapi_uri, ntemplates, structure_dirs, loopmodel):
    gather_templates_from_targetexplorer_metadata = gen_metadata_gather_from_targetexplorer(search_string, dbapi_uri)
    gather_templates_from_targetexplorer_metadata['structure_dirs'] = structure_dirs
    gather_templates_metadata = gen_gather_templates_metadata(ntemplates, loopmodel=loopmodel, additional_metadata=gather_templates_from_targetexplorer_metadata)
    project_metadata = msmseeder.core.ProjectMetadata(project_stage='gather_templates')
    project_metadata.add_data(gather_templates_metadata)
    project_metadata.write()


@msmseeder.utils.mpirank0only_and_end_with_barrier
def write_gather_templates_from_uniprot_metadata(uniprot_query_string, uniprot_domain_regex, ntemplates, structure_dirs, loopmodel):
    gather_templates_from_uniprot_metadata = gen_metadata_gather_from_uniprot(uniprot_query_string, uniprot_domain_regex)
    gather_templates_from_uniprot_metadata['structure_dirs'] = structure_dirs
    gather_templates_metadata = gen_gather_templates_metadata(ntemplates, loopmodel=loopmodel, additional_metadata=gather_templates_from_uniprot_metadata)
    project_metadata = msmseeder.core.ProjectMetadata(project_stage='gather_templates')
    project_metadata.add_data(gather_templates_metadata)
    project_metadata.write()


def gen_gather_templates_metadata(nselected_templates, loopmodel=False, additional_metadata=None):
    if additional_metadata is None:
        additional_metadata = {}
    datestamp = msmseeder.core.get_utcnow_formatted()
    metadata = {
        'datestamp': datestamp,
        'ntemplates': str(nselected_templates),
        'loopmodel': loopmodel,
        'python_version': sys.version.split('|')[0].strip(),
        'python_full_version': msmseeder.core.literal_str(sys.version),
        'msmseeder_version': msmseeder.version.short_version,
        'msmseeder_commit': msmseeder.version.git_revision
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
                extensions={(None, 'match_regex'): msmseeder.core.xpath_match_regex_case_sensitive}
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
                    chain_spans = msmseeder.UniProt.parse_uniprot_pdbref_chains(chain_span_string)

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


def pdbfix_templates(selected_templates):
    missing_residues_sublist = []
    ntemplates = len(selected_templates)
    for template_index in range(mpistate.rank, ntemplates, mpistate.size):
        print template_index
        missing_residues_sublist.append(pdbfix_template(selected_templates[template_index]))

    missing_residues_gathered = mpistate.comm.gather(missing_residues_sublist, root=0)

    missing_residues_list = []
    if mpistate.rank == 0:
        missing_residues_list = [None] * ntemplates
        for template_index in range(ntemplates):
            missing_residues_list[template_index] = missing_residues_gathered[template_index % mpistate.size][template_index // mpistate.size]

    missing_residues_list = mpistate.comm.bcast(missing_residues_list, root=0)

    return missing_residues_list


def pdbfix_template(template):
    try:
        import pdbfixer
        import simtk.openmm.app
        import Bio.SeqUtils
        template_filepath = os.path.join(msmseeder.core.default_project_dirnames.templates_structures_resolved, template.templateid + '.pdb')
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
        template_pdbfixed_filepath = os.path.join(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-pdbfixed.pdb')
        with open(template_pdbfixed_filepath, 'w') as template_pdbfixed_file:
            simtk.openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, file=template_pdbfixed_file)
        return fixer.missingResidues
    except KeyboardInterrupt:
        raise
    except Exception as e:
        trbk = traceback.format_exc()
        log_filepath = os.path.abspath(os.path.join(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-pdbfixer-log.yaml'))
        logfile = msmseeder.core.LogFile(log_filepath)
        logfile.log({
            'templateid': str(template.templateid),
            'exception': e,   # if Rosetta loopmodel fails, this field should include text sent to stdout and stderr
            'traceback': msmseeder.core.literal_str(trbk),
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


def loopmodel_templates(selected_templates, missing_residues):
    for template_index in range(mpistate.rank, len(selected_templates), mpistate.size):
        template = selected_templates[template_index]
        if mpistate.size > 1:
            logger.info('MPI rank %d modeling missing loops for template %s' % (mpistate.rank, template.templateid))
        else:
            logger.info('Modeling missing loops for template %s' % template.templateid)
        loopmodel_template(template, missing_residues[template_index])


def loopmodel_template(template, missing_residues):
    template_filepath = os.path.abspath(os.path.join(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-pdbfixed.pdb'))
    output_pdb_filepath = os.path.abspath(os.path.join(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '.pdb'))
    loop_filepath = os.path.abspath(os.path.join(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '.loop'))
    output_score_filepath = os.path.abspath(os.path.join(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-loopmodel-score.sc'))
    log_filepath = os.path.abspath(os.path.join(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '-loopmodel-log.yaml'))
    logfile = msmseeder.core.LogFile(log_filepath)
    try:
        write_loop_file(template, missing_residues)
        if len(missing_residues) == 0:
            shutil.copy(template_filepath, output_pdb_filepath)
            return
        loopmodel_output_text = run_loopmodel(template_filepath, loop_filepath, output_pdb_filepath, output_score_filepath)
        logfile.log({
            'templateid': str(template.templateid),
            'loopmodel_output': msmseeder.core.literal_str(loopmodel_output_text),
            'mpi_rank': mpistate.rank,
            'successful': True,
            })
    except KeyboardInterrupt:
        raise
    except Exception as e:
        trbk = traceback.format_exc()
        logfile.log({
            'templateid': str(template.templateid),
            'exception': e,   # if Rosetta loopmodel fails, this field should include text sent to stdout and stderr
            'traceback': msmseeder.core.literal_str(trbk),
            'mpi_rank': mpistate.rank,
            'successful': False,
        })
        logger.error('MPI rank %d Loopmodel error for template %s - see logfile' % (mpistate.rank, template.templateid))


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
    loop_filepath = os.path.join(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops, template.templateid + '.loop')
    with open(loop_filepath, 'w') as loop_file:
        loop_file.write(loop_file_text)


def run_loopmodel(input_template_pdb_filepath, loop_filepath, output_pdb_filepath, output_score_filepath, loopmodel_executable_filepath=msmseeder.core.find_loopmodel_executable(), nmodels_to_build=1):
    temp_dir = tempfile.mkdtemp()
    temp_template_filepath = os.path.join(temp_dir, 'template.pdb')
    temp_loop_filepath = os.path.join(temp_dir, 'template.loop')
    try:
        minirosetta_database_path = os.environ.get('MINIROSETTA_DATABASE')
        shutil.copy(input_template_pdb_filepath, temp_template_filepath)
        shutil.copy(loop_filepath, temp_loop_filepath)
        output_text = subprocess.check_output(
            [
                loopmodel_executable_filepath,
                '-database', minirosetta_database_path,
                '-in::file::s', temp_template_filepath,
                '-loops:loop_file', temp_loop_filepath,
                '-loops:remodel', 'perturb_kic',
                '-loops:refine', 'refine_kic',
                '-ex1',
                '-ex2',
                '-nstruct', '%d' % nmodels_to_build,
                '-loops:max_kic_build_attempts', '100',
                '-in:file:fullatom',
            ],
            stderr=subprocess.STDOUT
        )
        shutil.copy('template_0001.pdb', output_pdb_filepath)
        shutil.copy('score.sc', output_score_filepath)
        shutil.rmtree(temp_dir)
        return output_text
    except:
        # os.chdir(cwd)
        shutil.rmtree(temp_dir)
        raise














        # def create_dummy_complete_templates(selected_templates):
        #     # TODO
        #     for template in selected_templates:
        #         template_pdbresnums = template['template_pdbresnums']
        #         template_seq_complete = template['template_seq']


        # def pdbfix_template(template):
        #     # TODO
        #     import pdbfixer
        #     import simtk.openmm.app as app
        #
        #     missing_residues = determine_template_missing_residues(template.resolved_seq, template.full_seq)
        #
        #     pdbfixer.findMissingAtoms()
        #     pdbfixer.addMissingAtoms()
        #     app.PDBFile.writeFile(pdbfixer.topology, pdbfixer.positions, file=template_fixed_filepath)


        # def determine_template_missing_residues(template):
        #     # full_seq_residue_indices = range(len(template.resolved_seq))
        #
        #     # missing_residues = {}
        #     # for i in range(len(template.full_seq)):
        #     #     if template.full_pdbresnums[i] not in template.resolved_pdbresnums[i]:
        #     #         missing_residues[(0, i)] = []
        #
        #     # (0,17): [
        #     #         'GLY',
        #     #         'SER',
        #     #         'PHE',
        #     #         'GLY',
        #     #      ]
        #
        #     # seq_observed = template.resolved_seq
        #     seq_complete = template.full_seq
        #     seq_observed_w_gaps = gen_seq_observed_w_gaps(template)
        #     seq_observed_offset = index of first non-null element in seq_observed_w_gaps
        #
        #     missing_residues = {}
        #     index = 0
        #     for i in range(len(seq_complete)):
        #         if i < seq_observed_offset or i >= len(seq_observed_w_gaps)+seq_observed_offset or seq_observed_w_gaps[i-seq_observed_offset] is None:
        #             key = (0, index)
        #             if key not in missing_residues:
        #                 missing_residues[key] = []
        #             residue_code = seq_complete[i]
        #             residue_name = Bio.SeqUtils.seq3(residue_code).upper()
        #             missing_residues[key].append(residue_name)
        #         else:
        #             index += 1


        # def gen_seq_observed_w_gaps(template):
        #     # TODO skip this - should be able to set the "complete" sequence with the PDBFixer object and then just use the findMissingResidues() function.
        #     first_observed_residue_index = template.full_pdbresnums.index(template.resolved_pdbresnums[0])
        #     last_observed_residue_index = template.full_pdbresnums.index(template.resolved_pdbresnums[-1])
        #     nobserved_residues_w_gaps = last_observed_residue_index - first_observed_residue_index + 1
        #     print first_observed_residue_index, last_observed_residue_index, nobserved_residues_w_gaps, len(template.full_seq)
        #
        #     seq_observed_w_gaps = [None] * nobserved_residues_w_gaps
        #     for i in range(len(template.full_seq)):
        #         if i >= first_observed_residue_index and i <= last_observed_residue_index:
        #             index = i - first_observed_residue_index
        #             seq_observed_w_gaps[index] = template.full_seq[index]
        #     return seq_observed_w_gaps