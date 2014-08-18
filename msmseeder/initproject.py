import numpy as np

# =========
# Project initialization
# =========

def init(project_toplevel_dir):
    '''Initialize MSMSeeder project within the current directory. Creates
    necessary subdirectories and a project metadata .yaml file.
    '''

    # =========
    # Parameters
    # =========

    import sys
    import os
    import datetime
    import json
    import yaml
    import msmseeder
    import msmseeder.version

    project_dirnames = ['targets', 'structures', 'templates', 'models', 'packaged-models']

    datestamp = msmseeder.core.get_utcnow_formatted()

    # =========
    # Create necessary project directories
    # =========

    os.chdir(project_toplevel_dir)

    for dirname in project_dirnames:
        try:
            os.mkdir(dirname)
            print 'Created directory "%s"' % dirname
        except OSError as e:
            if e.errno == 17:
                print 'Directory "%s" already exists - will not overwrite' % dirname
            else:
                raise
    os.mkdir(os.path.join('structures', 'pdb'))
    os.mkdir(os.path.join('structures', 'sifts'))
    os.mkdir(os.path.join('templates', 'structures'))

    # =========
    # Metadata
    # =========

    metadata = {
        'init': {
            'datestamp': datestamp,
            'init_path': os.path.abspath(project_toplevel_dir),
            'python_version': sys.version.split('|')[0].strip(),
            'python_full_version': sys.version,
            'msmseeder_version': msmseeder.version.short_version,
            'msmseeder_commit': msmseeder.version.git_revision
        }
    }

    metadata = msmseeder.core.ProjectMetadata(metadata)
    metadata.write('meta.yaml')

    print 'Done.'


# =========
# Gather targets methods
# =========

def gather_targets_from_TargetExplorerDB(DB_path, species=None):
    '''Gather protein target data from an existing TargetExplorerDB database.'''

    # =========
    # Parameters
    # =========

    import sys
    import os
    import datetime
    import imp
    import yaml
    import msmseeder
    import msmseeder.version
    from lxml import etree

    fasta_ofilepath = os.path.join('targets', 'targets.fa')

    # =========
    # Parse the TargetExplorer database path
    # =========

    parser = etree.XMLParser(huge_tree=True)
    DB_root = etree.parse(DB_path, parser).getroot()

    # =========
    # Extract target data from database
    # =========

    print 'Extracting target data from database...'

    if species != None:
        target_domains = DB_root.findall('entry/UniProt[@NCBI_taxID="%s"]/domains/domain[@targetID]' % species)
    else:
        target_domains = DB_root.findall('entry/UniProt/domains/domain[@targetID]')

    # =========
    # Also get the original uniprot search strings from the database project directory
    # =========

    db_config_path = os.path.abspath(os.path.join(os.path.dirname(DB_path), '..', 'project_config.py'))
    target_explorer_db_project_config = imp.load_source('project_config', db_config_path)

    # =========
    # Write target data to FASTA file
    # =========

    print 'Writing target data to FASTA file "%s"...' % fasta_ofilepath

    with open(fasta_ofilepath, 'w') as fasta_ofile:
        for target_domain in target_domains:
            targetID = target_domain.get('targetID')
            targetseqnode = target_domain.find('sequence')
            targetseq = targetseqnode.text.strip()

            target_fasta_string = '>%s\n%s\n' % (targetID, targetseq)

            fasta_ofile.write(target_fasta_string)

    # =========
    # Metadata
    # =========

    datestamp = msmseeder.core.get_utcnow_formatted()

    with open('meta.yaml') as init_meta_file:
        metadata = yaml.load(init_meta_file)

    metadata['gather_targets'] = {
        'datestamp': datestamp,
        'method': 'TargetExplorerDB',
        'ntargets': str(len(target_domains)),
        'gather_from_target_explorer': {
            'uniprot_query_string': target_explorer_db_project_config.uniprot_query_string,
            'uniprot_domain_regex': target_explorer_db_project_config.uniprot_domain_regex,
            'species_selector': species if not None else '',
            'database_path': DB_path
        },
        'python_version': sys.version.split('|')[0].strip(),
        'python_full_version': sys.version,
        'msmseeder_version': msmseeder.version.short_version,
        'msmseeder_commit': msmseeder.version.git_revision
    }

    metadata = msmseeder.core.ProjectMetadata(metadata)
    metadata.write('targets/meta.yaml')

    print 'Done.'

def gather_targets_from_UniProt():
    '''# Searches UniProt for a set of target proteins with a user-defined
    query string, then saves target IDs and sequences.'''

    raise Exception, 'Not yet implemented.'

    # TODO NOTE use 'ABL1_HUMAN_D0' style for targetIDs in this script. Don't bother with mutants. The gather_targets_from_TargetExplorerDB version of this routine will use whatever targetIDs it finds in the DB, which may include 'ABL1_HUMAN_D0_M0' for example.
    # TODO NOTE read in manual exceptions file. Should at least handle manual exceptions for domain spans (e.g. for ABL1_HUMAN_D0)

    # =========
    # Parameters
    # =========

    # =========
    # Get user-defined UniProt query string
    # =========

    # check for command-line arg

    # and check the project metadata file


# =========
# Gather templates methods
# =========

def gather_templates_from_TargetExplorerDB(DB_path):
    '''Gather protein target data from an existing TargetExplorerDB database.'''

    raise Exception, 'Not implemented yet.'

    # =========
    # Parameters
    # =========

    import os, datetime
    import msmseeder
    from lxml import etree

    fasta_ofilepath = os.path.join('templates', 'templates.fa')

    # =========
    # Read in project metadata file
    # =========

    project_metadata = msmseeder.core.ProjectMetadata()
    project_metadata.load(msmseeder.core.project_metadata_filename)

    # =========
    # Parse the TargetExplorer database
    # =========

    # Parse the DB
    parser = etree.XMLParser(huge_tree=True)
    DB_root = etree.parse(DB_path, parser).getroot()

    # =========
    # Extract template IDs data from database
    # =========

    print 'Extracting template data from database...'

    # TODO

    # =========
    # Download structures if necessary
    # =========

    # TODO

    # =========
    # Write template data to FASTA file
    # =========

    print 'Writing template data to FASTA file "%s"...' % fasta_ofilepath

    # TODO

    # with open(fasta_ofilepath, 'w') as fasta_ofile:
    #     for template in templates:
    #         templateID =
    #         templateseq =
    #
    #         template_fasta_string = '>%s\n%s\n' % (templateID, templateseq)
    #
    #         fasta_ofile.write(template_fasta_string)

    # =========
    # Update project metadata file
    # =========

    datestamp = msmseeder.core.get_utcnow_formatted()

    template_selection_metadata = {
    'datestamp': datestamp,
    'template-selection-method': 'TargetExplorerDB',
    'TargetExplorerDB-database-path': DB_path
    }

    project_metadata.data['template-selection'] = template_selection_metadata
    project_metadata.write()

    print 'Done.'


def gather_templates_from_UniProt(UniProt_query_string, UniProt_domain_regex, structure_paths=[]):
    '''# Searches UniProt for a set of template proteins with a user-defined
    query string, then saves IDs, sequences and structures.'''

    # =========
    # Parameters
    # =========

    import sys
    import os
    import datetime
    import yaml
    import gzip
    import msmseeder
    import msmseeder.UniProt
    import msmseeder.PDB
    import msmseeder.version
    from lxml import etree

    fasta_ofilepath = os.path.join('templates', 'templates.fa')

    # =========
    # Read in project manual specifications
    # =========

    min_domain_len = None
    max_domain_len = None
    domain_span_manual_specifications = None
    skip_pdbs = []
    if os.path.exists(msmseeder.core.manual_specifications_filename):
        with open(msmseeder.core.manual_specifications_filename, 'r') as manual_specifications_file:
            manual_specifications = yaml.load(manual_specifications_file)
        template_manual_specifications = manual_specifications.get('template-selection')
        min_domain_len = template_manual_specifications.get('min-domain-len') if not None else 0
        max_domain_len = template_manual_specifications.get('max-domain-len') if not None else np.nan
        domain_span_manual_specifications = template_manual_specifications.get('domain-spans') if not None else {}
        skip_pdbs = template_manual_specifications.get('skip-pdbs') if not None else []

    # =========
    # Make request to UniProt web server and parse the returned XML
    # =========

    print 'Querying UniProt web server...'
    UniProtXMLstring = msmseeder.UniProt.retrieve_uniprot(UniProt_query_string)
    parser = etree.XMLParser(huge_tree=True)
    UniProtXML = etree.fromstring(UniProtXMLstring, parser)
    print 'Number of entries returned from initial UniProt search:', len(UniProtXML)
    print ''

    # =========
    # If the UniProt query string contained a domain selector, print the set of
    # unique UniProt domain names which would have been selected during the
    # UniProt search (case-insensitive). This should aid users in the
    # construction an appropriate regex for selecting an appropriate subset of
    # these domains.
    # =========

    if 'domain:' in UniProt_query_string:

        # First extract the domain selection
        # Example query string: 'domain:"Protein kinase" AND reviewed:yes'
        # Can assume that the domain selection will be bounded by double-quotes
        query_string_split = UniProt_query_string.split('"')
        query_string_domain_selection = query_string_split[ query_string_split.index('domain:') + 1 ]

        UniProt_query_string_domains = UniProtXML.xpath('entry/feature[@type="domain"][match_regex(@description, "%s")]' % query_string_domain_selection, extensions = { (None, 'match_regex'): msmseeder.core.xpath_match_regex_case_insensitive })

        UniProt_unique_domain_names = set([domain.get('description') for domain in UniProt_query_string_domains])
        print 'Set of unique domain names selected by the domain selector \'%s\' during the initial UniProt search:\n%s' % (query_string_domain_selection, UniProt_unique_domain_names)
        print ''

    else:
        UniProt_domains = UniProtXML.xpath('entry/feature[@type="domain"]')
        UniProt_unique_domain_names = set([domain.get('description') for domain in UniProt_domains])
        print 'Set of unique domain names returned from the initial UniProt search using the query string \'%s\':\n%s' % (UniProt_query_string, UniProt_unique_domain_names)
        print ''

    # =========
    # Print subset of domains returned following filtering with the UniProt_domain_regex (case sensitive)
    # =========

    if UniProt_domain_regex != None:
        regex_matched_domains = UniProtXML.xpath('entry/feature[@type="domain"][match_regex(@description, "%s")]' % UniProt_domain_regex, extensions = { (None, 'match_regex'): msmseeder.core.xpath_match_regex_case_sensitive })

        regex_matched_domains_unique_names = set([domain.get('description') for domain in regex_matched_domains])
        print 'Unique domain names selected after searching with the case-sensitive regex string \'%s\':\n%s' % (UniProt_domain_regex, regex_matched_domains_unique_names)
        print ''

    # =========
    # Now go through all the UniProt entries, domains, PDBs and chains, do some filtering, and store data for the selected PDB chains
    # =========

    all_UniProt_entries = UniProtXML.findall('entry')

    selected_PDBchains = []

    print 'Extracting information from returned UniProt data...'

    for entry in all_UniProt_entries:
        entry_name = entry.find('name').text
        if UniProt_domain_regex != None:
            selected_domains = entry.xpath('feature[@type="domain"][match_regex(@description, "%s")]' % UniProt_domain_regex, extensions = { (None, 'match_regex'): msmseeder.core.xpath_match_regex_case_sensitive })
        else:
            selected_domains = entry.findall('feature[@type="domain"]')

        domain_iter = 0
        for domain in selected_domains:
            domainID = '%(entry_name)s_D%(domain_iter)d' % vars()
            domain_span = [int(domain.find('location/begin').get('position')), int(domain.find('location/end').get('position'))]
            if domainID in domain_span_manual_specifications:
                domain_span = [int(x) for x in domain_span_manual_specifications[domainID].split('-')]
            domain_len = domain_span[1] - domain_span[0] + 1
            if domain_len < min_domain_len or domain_len > max_domain_len:
                continue

            domain_iter += 1
            PDBs = domain.getparent().xpath('dbReference[@type="PDB"]/property[@type="method"][@value="X-ray" or @value="NMR"]/..')

            for PDB in PDBs:
                PDBID = PDB.get('id')
                if PDBID in skip_pdbs:
                    continue
                PDB_chain_span_nodes = PDB.findall('property[@type="chains"]')

                for PDB_chain_span_node in PDB_chain_span_nodes:
                    chain_span_string = PDB_chain_span_node.get('value')
                    chain_spans = msmseeder.UniProt.parse_uniprot_pdbref_chains(chain_span_string)

                    for chainID in chain_spans.keys():
                        span = chain_spans[chainID]
                        if (span[0] < domain_span[0]+30) & (span[1] > domain_span[1]-30):
                            templateID = '%(domainID)s_%(PDBID)s_%(chainID)s' % vars()
                            data = {
                            'templateID': templateID,
                            'PDBID': PDBID,
                            'chainID': chainID,
                            'domain_span': domain_span
                            }
                            selected_PDBchains.append(data)

    print '%d PDB chains selected.' % len(selected_PDBchains)
    print ''

    # =========
    # Search for PDB and SIFTS files; download if necessary
    # =========

    for PDBchain in selected_PDBchains:
        PDBID = PDBchain['PDBID']
        chainID = PDBchain['chainID']

        project_pdb_filepath = os.path.join('structures', 'pdb', PDBID + '.pdb')
        project_sifts_filepath = os.path.join('structures', 'sifts', PDBID + '.xml.gz')

        # Check if PDB file/symlink already exists and is not empty
        search_for_pdb = True
        if os.path.exists(project_pdb_filepath):
            if os.path.getsize(project_pdb_filepath) > 0:
                search_for_pdb = False

        # If not, search any user-defined paths and create a symlink if found
        if search_for_pdb:
            for structure_dir in structure_paths:
                pdb_filepath = os.path.join(structure_dir, PDBID + '.pdb')
                if os.path.exists(pdb_filepath):
                    if os.path.getsize(pdb_filepath) > 0:
                        if os.path.exists(project_pdb_filepath):
                            os.remove(project_pdb_filepath)
                        os.symlink(pdb_filepath, project_pdb_filepath)
                        break

            # If still not found, download the PDB file
            if not os.path.exists(project_pdb_filepath):
                print 'Downloading PDB file for:', PDBID
                pdbgz_page = msmseeder.PDB.retrieve_pdb(PDBID, compressed='yes')
                with open(project_pdb_filepath + '.gz', 'w') as pdbgz_file:
                    pdbgz_file.write(pdbgz_page)    
                with gzip.open(project_pdb_filepath + '.gz', 'rb') as pdbgz_file_decoded:
                    with open(project_pdb_filepath, 'w') as project_pdb_file:
                        project_pdb_file.writelines(pdbgz_file_decoded)
                os.remove(project_pdb_filepath + '.gz')

        # Check if SIFTS file already exists and is not empty
        search_for_sifts = True
        if os.path.exists(project_sifts_filepath):
            if os.path.getsize(project_sifts_filepath) > 0:
                search_for_sifts = False

        # If not, search any user-defined paths and create a symlink if found
        if search_for_sifts:
            for structure_dir in structure_paths:
                sifts_filepath = os.path.join(structure_dir, PDBID + '.xml.gz')
                if os.path.exists(sifts_filepath):
                    if os.path.getsize(sifts_filepath) > 0:
                        if os.path.exists(project_sifts_filepath):
                            os.remove(project_sifts_filepath)
                        os.symlink(sifts_filepath, project_sifts_filepath)
                        break

            # If still not found, download the SIFTS file
            if not os.path.exists(project_sifts_filepath):
                print 'Downloading sifts file for:', PDBID
                sifts_page = msmseeder.PDB.retrieve_sifts(PDBID)
                with gzip.open(project_sifts_filepath, 'wb') as project_sifts_file: 
                    project_sifts_file.write(sifts_page)    

    # =========
    # Extract PDBchain residues using SIFTS files
    # =========

    selected_templates = []

    print 'Extracting residues from PDB chains...'

    for PDBchain in selected_PDBchains:
        templateID = PDBchain['templateID']
        chainID = PDBchain['chainID']
        PDBID = PDBchain['PDBID']
        domain_span = PDBchain['domain_span']

        # parse SIFTS XML document
        sifts_filepath = os.path.join('structures', 'sifts', PDBID + '.xml.gz')
        with gzip.open(sifts_filepath, 'rb') as sifts_file:
            parser = etree.XMLParser(huge_tree=True)
            siftsXML = etree.parse(sifts_file, parser).getroot()

        # firstly, add "PDB modified" tags to certain phosphorylated residue types, which sometimes do not have such tags in the SIFTS file
        # known cases: 4BCP, 4BCG, 4I5C, 4IVB, 4IAC
        modified_residues = []
        modified_residues += siftsXML.findall('entity/segment/listResidue/residue[@dbResName="TPO"]')
        modified_residues += siftsXML.findall('entity/segment/listResidue/residue[@dbResName="PTR"]')
        modified_residues += siftsXML.findall('entity/segment/listResidue/residue[@dbResName="SEP"]')
        for mr in modified_residues:
            if mr == None:
                continue
            residue_detail_modified = etree.Element('residueDetail')
            residue_detail_modified.set('dbSource','MSD')
            residue_detail_modified.set('property','Annotation')
            residue_detail_modified.text = 'PDB\n          modified'
            mr.append(residue_detail_modified)

        # now extract PDB residues with the correct PDB chain ID, are observed, have a UniProt crossref and are within the UniProt domain bounds, and do not have "PDB modified", "Conflict" or "Engineered mutation" tags.
        selected_residues = siftsXML.xpath('entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][not(../residueDetail[contains(text(),"Not_Observed")])][../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]][not(../residueDetail[contains(text(),"modified")])][not(../residueDetail[contains(text(),"Conflict")])][not(../residueDetail[contains(text(),"mutation")])]' % (chainID, domain_span[0], domain_span[1]))
        # calculate the ratio of observed residues - if less than a certain amount, discard PDBchain
        all_PDB_domain_residues = siftsXML.xpath('entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][../crossRefDb[@dbSource="UniProt"][@dbResNum >= "%d"][@dbResNum <= "%d"]]' % (chainID, domain_span[0], domain_span[1]))
        if len(selected_residues) == 0 or len(all_PDB_domain_residues) == 0:
            continue

        ratio_observed = float(len(selected_residues)) / float(len(all_PDB_domain_residues))
        if ratio_observed < msmseeder.core.template_acceptable_ratio_observed_residues:
            #PDBchain['DISCARD'] = True
            continue

        # make a single-letter aa code sequence
        template_seq = ''.join([residue.find('../crossRefDb[@dbSource="UniProt"]').get('dbResName') for residue in selected_residues])
        # store as strings to match resnums in PDB file directly. Important because PDB resnums may be e.g. '56A' in 3HLL
        template_PDBresnums = [residue.get('dbResNum') for residue in selected_residues]

        # store data
        template_data = {
        'PDBID': PDBID,
        'chainID': chainID,
        'templateID': templateID,
        'template_seq': template_seq,
        'template_PDBresnums': template_PDBresnums
        }
        selected_templates.append(template_data)

    print '%d templates selected.' % len(selected_templates)
    print ''

    # =========
    # Write template IDs and sequences to file
    # =========

    print 'Writing template IDs and sequences to file:', fasta_ofilepath

    with open(fasta_ofilepath, 'w') as fasta_ofile:
        for template in selected_templates:
            fasta_ofile.write('>' + template['templateID'] + '\n')
            fasta_ofile.write(template['template_seq'] + '\n')

    # =========
    # Extract template structures from PDB files and write to file
    # =========

    print 'Writing template structures...'

    for template in selected_templates:
        PDBID = template['PDBID']
        chainID = template['chainID']
        templateID = template['templateID']
        template_PDBresnums = template['template_PDBresnums']
        pdb_filename = os.path.join('structures', 'pdb', PDBID + '.pdb')
        template_filename = os.path.join('templates', 'structures', templateID + '.pdb')
        nresidues_extracted = msmseeder.PDB.extract_residues_by_resnum(template_filename, pdb_filename, template_PDBresnums, chainID)
        if nresidues_extracted != len(template_PDBresnums):
            raise Exception, 'Number of residues extracted from PDB file does not match desired number of residues.'

    # =========
    # Metadata
    # =========

    datestamp = msmseeder.core.get_utcnow_formatted()

    with open('meta.yaml') as init_meta_file:
        metadata = yaml.load(init_meta_file)

    metadata['gather_templates'] = {
        'datestamp': datestamp,
        'method': 'UniProt',
        'ntemplates': str(len(selected_templates)),
        'gather_from_uniprot': {
            'uniprot_query_string': UniProt_query_string,
            'uniprot_domain_regex': UniProt_domain_regex if not None else ''
        },
        'structure_paths': structure_paths,
        'python_version': sys.version.split('|')[0].strip(),
        'python_full_version': sys.version,
        'msmseeder_version': msmseeder.version.short_version,
        'msmseeder_commit': msmseeder.version.git_revision
    }

    metadata = msmseeder.core.ProjectMetadata(metadata)
    metadata.write('templates/meta.yaml')

    print 'Done.'
