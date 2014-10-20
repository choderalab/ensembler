import gzip
import json
import sys
import os
import msmseeder.version
from lxml import etree
import msmseeder
import msmseeder.TargetExplorer
import msmseeder.UniProt
import msmseeder.PDB
from msmseeder.core import construct_fasta_str
from msmseeder.utils import file_exists_and_not_empty


@msmseeder.utils.notify_when_done
def initproject(project_toplevel_dir):
    """Initialize MSMSeeder project within the given directory. Creates
    necessary subdirectories and a project metadata .yaml file.
    :param project_toplevel_dir: str
    """
    create_dirs(project_toplevel_dir)
    init_metadata = gen_init_metadata(project_toplevel_dir)
    msmseeder.core.write_metadata(init_metadata, msmseeder_stage='init')


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

    gather_from_targetexplorer_metadata = gen_metadata_gather_from_targetexplorer(search_string, dbapi_uri)
    gather_targets_metadata = gen_gather_targets_metadata(len(targets), additional_metadata=gather_from_targetexplorer_metadata)
    msmseeder.core.write_metadata(gather_targets_metadata, msmseeder_stage='gather_targets')


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

    gather_from_uniprot_metadata = gen_metadata_gather_from_uniprot(uniprot_query_string, uniprot_domain_regex)
    gather_targets_metadata = gen_gather_targets_metadata(len(targets), additional_metadata=gather_from_uniprot_metadata)
    msmseeder.core.write_metadata(gather_targets_metadata, msmseeder_stage='gather_targets')


def create_dirs(project_toplevel_dir):
    os.chdir(project_toplevel_dir)
    for dirname in msmseeder.core.project_dirnames:
        msmseeder.utils.create_dir(dirname)


def gen_init_metadata(project_toplevel_dir):
    datestamp = msmseeder.core.get_utcnow_formatted()
    metadata_dict = {
        'init': {
            'datestamp': datestamp,
            'init_path': os.path.abspath(project_toplevel_dir),
            'python_version': sys.version.split('|')[0].strip(),
            'python_full_version': msmseeder.core.literal_str(sys.version),
            'msmseeder_version': msmseeder.version.short_version,
            'msmseeder_commit': msmseeder.version.git_revision
        }
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
    print 'Querying UniProt web server...'
    uniprotxmlstring = msmseeder.UniProt.retrieve_uniprot(uniprot_query_string)
    parser = etree.XMLParser(huge_tree=True)
    uniprotxml = etree.fromstring(uniprotxmlstring, parser)
    print 'Number of entries returned from initial UniProt search: %r\n' % len(uniprotxml)
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
    print 'Writing target data to FASTA file "%s"...' % fasta_ofilepath
    with open(fasta_ofilepath, 'w') as fasta_ofile:
        for target in targets:
            targetid, targetseq = target
            target_fasta_string = construct_fasta_str(targetid, targetseq)
            fasta_ofile.write(target_fasta_string)


def gen_metadata_gather_from_targetexplorer(search_string, dbapi_uri):
    db_metadata = get_db_metadata(dbapi_uri)
    metadata = {
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
        'gather_from_uniprot': {
            'uniprot_query_string': uniprot_query_string,
            'uniprot_domain_regex': uniprot_domain_regex,
        }
    }
    return metadata


def get_db_metadata(dbapi_uri):
    db_metadata_jsonstr = msmseeder.TargetExplorer.get_targetexplorer_metadata(dbapi_uri)
    return json.loads(db_metadata_jsonstr)


def gen_gather_targets_metadata(ntarget_domains, additional_metadata=None):
    if additional_metadata is None:
        additional_metadata = {}
    datestamp = msmseeder.core.get_utcnow_formatted()
    metadata = {
        'gather_targets': {
            'datestamp': datestamp,
            'method': 'TargetExplorerDB',
            'ntargets': str(ntarget_domains),
            'python_version': sys.version.split('|')[0].strip(),
            'python_full_version': msmseeder.core.literal_str(sys.version),
            'msmseeder_version': msmseeder.version.short_version,
            'msmseeder_commit': msmseeder.version.git_revision,
        }
    }
    metadata['gather_targets'].update(additional_metadata)
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
                (None, 'match_regex'): msmseeder.core.xpath_match_regex_case_insensitive
            }
        )
        uniprot_unique_domain_names = set([domain.get('description') for domain in uniprot_query_string_domains])
        print 'Set of unique domain names selected by the domain selector \'%s\'during the initial UniProt search:\n%s\n' \
              % (query_string_domain_selection, uniprot_unique_domain_names)

    else:
        uniprot_domains = uniprotxml.xpath('entry/feature[@type="domain"]')
        uniprot_unique_domain_names = set([domain.get('description') for domain in uniprot_domains])
        print 'Set of unique domain names returned from the initial UniProt search using the query string \'%s\':\n%s\n' \
              % (uniprot_query_string, uniprot_unique_domain_names)


def log_unique_domain_names_selected_by_regex(uniprot_domain_regex, uniprotxml):
    regex_matched_domains = uniprotxml.xpath(
        'entry/feature[@type="domain"][match_regex(@description, "%s")]' % uniprot_domain_regex,
        extensions={(None, 'match_regex'): msmseeder.core.xpath_match_regex_case_sensitive}
    )
    regex_matched_domains_unique_names = set([domain.get('description') for domain in regex_matched_domains])
    print 'Unique domain names selected after searching with the case-sensitive regex string \'%s\':\n%s\n' \
          % (uniprot_domain_regex, regex_matched_domains_unique_names)


@msmseeder.utils.notify_when_done
def gather_templates_from_targetexplorerdb(dbapi_uri, search_string='', structure_dirs=None):
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

    additional_metadata = gen_metadata_gather_from_targetexplorer(search_string, dbapi_uri)
    additional_metadata['structure_dirs'] = structure_dirs
    gather_templates_metadata = gen_gather_templates_metadata(len(selected_templates), additional_metadata)
    msmseeder.core.write_metadata(gather_templates_metadata, msmseeder_stage='gather_templates')


def get_targetexplorer_templates_json(dbapi_uri, search_string):
    """
    :param dbapi_uri: str
    :param search_string: str
    :return: list containing nested lists and dicts
    """
    targetexplorer_jsonstr = msmseeder.TargetExplorer.query_targetexplorer(
        dbapi_uri, search_string, return_data='pdb_data'
    )
    targetexplorer_json = json.loads(targetexplorer_jsonstr)
    return targetexplorer_json


def extract_template_pdbchains_from_targetexplorer_json(targetexplorer_json, manual_overrides):
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


def download_pdb_file(pdbid, project_pdb_filepath):
    print 'Downloading PDB file for:', pdbid
    pdbgz_page = msmseeder.PDB.retrieve_pdb(pdbid, compressed='yes')
    with open(project_pdb_filepath, 'w') as pdbgz_file:
        pdbgz_file.write(pdbgz_page)


def download_sifts_file(pdbid, project_sifts_filepath):
    print 'Downloading sifts file for:', pdbid
    sifts_page = msmseeder.PDB.retrieve_sifts(pdbid)
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
                structure_downloader = structure_downloader_mapper[structure_type]
                structure_downloader(pdbid, project_structure_filepath)


def extract_template_pdb_chain_residues(selected_pdbchains):
    print 'Extracting residues from PDB chains...'
    selected_templates = []
    for pdbchain in selected_pdbchains:
        extracted_pdb_template_seq_data = extract_pdb_template_seq(pdbchain)
        if extracted_pdb_template_seq_data is not None:
            selected_templates.append(extracted_pdb_template_seq_data)
    print '%d templates selected.\n' % len(selected_templates)
    return selected_templates


def extract_pdb_template_seq(pdbchain):
    """Extract data from PDB chain"""

    templateid = pdbchain['templateid']
    chainid = pdbchain['chainid']
    pdbid = pdbchain['pdbid']
    domain_span = pdbchain['domain_span']

    # parse SIFTS XML document
    sifts_filepath = os.path.join('structures', 'sifts', pdbid + '.xml.gz')
    with gzip.open(sifts_filepath, 'rb') as sifts_file:
        parser = etree.XMLParser(huge_tree=True)
        siftsxml = etree.parse(sifts_file, parser).getroot()

    # firstly, add "PDB modified" tags to certain phosphorylated residue types, which sometimes do not have such tags in the SIFTS file
    # known cases: 4BCP, 4BCG, 4I5C, 4IVB, 4IAC
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

    # now extract PDB residues with the correct PDB chain ID, are observed, have a UniProt crossref and are within the UniProt domain bounds, and do not have "PDB modified", "Conflict" or "Engineered mutation" tags.
    selected_residues = siftsxml.xpath(
        'entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][not(../residueDetail[contains(text(),"Not_Observed")])]'
        '[../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]][not(../residueDetail[contains(text(),"modified")])]'
        '[not(../residueDetail[contains(text(),"Conflict")])][not(../residueDetail[contains(text(),"mutation")])]' % (chainid, domain_span[0], domain_span[1])
    )
    # calculate the ratio of observed residues - if less than a certain amount, discard pdbchain
    all_pdb_domain_residues = siftsxml.xpath(
        'entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]]' % (chainid, domain_span[0], domain_span[1])
    )
    if len(selected_residues) == 0 or len(all_pdb_domain_residues) == 0:
        return

    ratio_observed = float(len(selected_residues)) / float(len(all_pdb_domain_residues))
    if ratio_observed < msmseeder.core.template_acceptable_ratio_observed_residues:
        return

    # make a single-letter aa code sequence
    template_seq = ''.join([residue.find('../crossRefDb[@dbSource="UniProt"]').get('dbResName') for residue in selected_residues])
    # store as strings to match resnums in PDB file directly. Important because PDB resnums may be e.g. '56A' in 3HLL
    template_pdbresnums = [residue.get('dbResNum') for residue in selected_residues]

    # store data
    template_data = {
    'pdbid': pdbid,
    'chainid': chainid,
    'templateid': templateid,
    'template_seq': template_seq,
    'template_pdbresnums': template_pdbresnums
    }

    return template_data


def write_template_seqs_to_fasta_file(selected_templates):
    selected_template_seq_tuples = [(template['templateid'], template['template_seq']) for template in selected_templates]
    write_seqs_to_fasta_file(selected_template_seq_tuples, fasta_ofilepath=os.path.join('templates', 'templates.fa'))


def extract_template_structures_from_pdb_files(selected_templates):
    print 'Writing template structures...'
    for template in selected_templates:
        pdbid = template['pdbid']
        chainid = template['chainid']
        templateid = template['templateid']
        template_pdbresnums = template['template_pdbresnums']
        pdb_filename = os.path.join('structures', 'pdb', pdbid + '.pdb.gz')
        template_filename = os.path.join('templates', 'structures', templateid + '.pdb')
        nresidues_extracted = msmseeder.PDB.extract_residues_by_resnum(
            template_filename, pdb_filename, template_pdbresnums, chainid
        )
        if nresidues_extracted != len(template_pdbresnums):
            raise Exception(
                'Number of residues extracted from PDB file (%d) does not match desired number of residues (%d).' % (
                    nresidues_extracted, template_pdbresnums
                )
            )


def gen_gather_templates_metadata(nselected_templates, additional_metadata=None):
    if additional_metadata is None:
        additional_metadata = {}
    datestamp = msmseeder.core.get_utcnow_formatted()
    metadata = {
        'gather_templates': {
            'datestamp': datestamp,
            'method': 'TargetExplorerDB',
            'ntemplates': str(nselected_templates),
            'python_version': sys.version.split('|')[0].strip(),
            'python_full_version': msmseeder.core.literal_str(sys.version),
            'msmseeder_version': msmseeder.version.short_version,
            'msmseeder_commit': msmseeder.version.git_revision
        }
    }
    metadata['gather_templates'].update(additional_metadata)
    return metadata


@msmseeder.utils.notify_when_done
def gather_templates_from_uniprot(uniprot_query_string, uniprot_domain_regex, structure_dirs=None):
    """# Searches UniProt for a set of template proteins with a user-defined
    query string, then saves IDs, sequences and structures."""
    manual_overrides = msmseeder.core.ManualOverrides()
    uniprotxml = get_uniprot_xml(uniprot_query_string)
    log_unique_domain_names(uniprot_query_string, uniprotxml)
    if uniprot_domain_regex is not None:
        log_unique_domain_names_selected_by_regex(uniprot_domain_regex, uniprotxml)

    selected_pdbchains = extract_template_pdbchains_from_uniprot_xml(uniprotxml, uniprot_domain_regex, manual_overrides)
    for pdbchain in selected_pdbchains:
        get_pdb_and_sifts_files(pdbchain['pdbid'], structure_dirs)

    selected_templates = extract_template_pdb_chain_residues(selected_pdbchains)
    write_template_seqs_to_fasta_file(selected_templates)
    extract_template_structures_from_pdb_files(selected_templates)

    additional_metadata = gen_metadata_gather_from_uniprot(uniprot_query_string, uniprot_domain_regex)
    additional_metadata['structure_dirs'] = structure_dirs
    gather_templates_metadata = gen_gather_templates_metadata(len(selected_templates), additional_metadata)
    msmseeder.core.write_metadata(gather_templates_metadata, msmseeder_stage='gather_templates')


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
    print '%d PDB chains selected.' % len(selected_pdbchains)
    print ''
    return selected_pdbchains


structure_type_file_extension_mapper = {'pdb': '.pdb.gz', 'sifts': '.xml.gz'}
structure_downloader_mapper = {'pdb': download_pdb_file, 'sifts': download_sifts_file}