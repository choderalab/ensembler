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

    import os, datetime, argparse
    import MSMSeeder

    project_dirnames = ['targets', 'templates', 'models', 'packaged-models']
    MSMSeeder.core.project_metadata_filename = 'project-data.yaml'

    now = datetime.datetime.utcnow()
    datestamp = now.strftime(MSMSeeder.core.datestamp_format_string)

    argparser = argparse.ArgumentParser(description='Initialize MSMSeeder project within the current directory. Creates necessary subdirectories and a project metadata .yaml file.')
    argparser.parse_args()

    # =========
    # Create necessary project directories
    # =========

    os.chdir(project_toplevel_dir)

    for dirname in project_dirnames:
        try:
            os.mkdir(dirname)
            print 'Created directory "%s"' % dirname
        except OSError as e:
            # If directory already exists, e.errno will be set to int 17, and e.strerror will be set to int 'File exists'
            if e.errno == 17:
                print 'Directory "%s" already exists - will not overwrite' % dirname
            else:
                raise
        
    # =========
    # Create metadata file and add datestamp
    # =========

    project_metadata = MSMSeeder.core.ProjectMetadata()
    init_metadata = {'init' : {'datestamp' : datestamp}}
    project_metadata.add_metadata(init_metadata)
    project_metadata.write(ofilepath=MSMSeeder.core.project_metadata_filename)

    print 'Done.'


# =========
# Gather targets methods
# =========

def gather_targets_from_TargetExplorerDB(DB_path):
    '''Gather protein target data from an existing TargetExplorerDB database.'''

    # =========
    # Parameters
    # =========

    import os, datetime
    import MSMSeeder
    from lxml import etree

    fasta_ofilepath = os.path.join('targets', 'targets.fa')

    # =========
    # Read in project metadata file (will throw an exception if it does not exist)
    # =========

    project_metadata = MSMSeeder.core.ProjectMetadata()
    project_metadata.load(MSMSeeder.core.project_metadata_filename)

    # =========
    # Parse the TargetExplorer database path
    # =========

    DB_root = etree.parse(DB_path).getroot()

    # =========
    # Extract target data from database
    # =========

    print 'Extracting target data from database...'

    target_domains = DB_root.findall('entry/UniProt/domains/domain[@targetID]')

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
    # Update project metadata file
    # =========

    now = datetime.datetime.utcnow()
    datestamp = now.strftime(MSMSeeder.core.datestamp_format_string)

    target_selection_metadata = {
    'datestamp': datestamp,
    'target-selection-method': 'TargetExplorerDB',
    'TargetExplorer-database-path': DB_path
    }

    project_metadata.data['target-selection'] = target_selection_metadata
    project_metadata.write()

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
    import MSMSeeder
    from lxml import etree

    fasta_ofilepath = os.path.join('templates', 'templates.fa')

    # =========
    # Read in project metadata file
    # =========

    project_metadata = MSMSeeder.core.ProjectMetadata()
    project_metadata.load(MSMSeeder.core.project_metadata_filename)

    # =========
    # Parse the TargetExplorer database
    # =========

    # Parse the DB
    DB_root = etree.parse(DB_path).getroot()
    DB_root

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

    now = datetime.datetime.utcnow()
    datestamp = now.strftime(MSMSeeder.core.datestamp_format_string)

    template_selection_metadata = {
    'datestamp': datestamp,
    'template-selection-method': 'TargetExplorerDB',
    'TargetExplorerDB-database-path': DB_path
    }

    project_metadata.data['template-selection'] = template_selection_metadata
    project_metadata.write()

    print 'Done.'


def gather_templates_from_UniProt(UniProt_query_string, UniProt_domain_regex):
    '''# Searches UniProt for a set of template proteins with a user-defined
    query string, then saves IDs, sequences and structures.'''

    # =========
    # Parameters
    # =========

    import os, datetime, yaml
    import MSMSeeder
    import MSMSeeder.UniProt
    from lxml import etree

    fasta_ofilepath = os.path.join('templates', 'templates.fa')

    # =========
    # Read in project metadata
    # =========

    project_metadata = MSMSeeder.core.ProjectMetadata()
    project_metadata.load(MSMSeeder.core.project_metadata_filename)

    if os.path.exists(MSMSeeder.core.manual_specifications_filename):
        with open(MSMSeeder.core.manual_specifications_filename, 'r') as manual_specifications_file:
            manual_specifications = yaml.load(manual_specifications_file)
        template_manual_specifications = manual_specifications.get('template-selection')
        min_domain_len = template_manual_specifications.get('min-domain-len')
        max_domain_len = template_manual_specifications.get('max-domain-len')
        domain_span_manual_specifications = template_manual_specifications.get('domain-spans')

    if min_domain_len == None:
        min_domain_len = 0
    if max_domain_len == None:
        max_domain_len = np.nan
    if domain_span_manual_specifications == None:
        domain_span_manual_specifications = {}

    # =========
    # Make request to UniProt web server and parse the returned XML
    # =========

    print 'Querying UniProt web server...'
    UniProtXMLstring = MSMSeeder.UniProt.retrieve_uniprot(UniProt_query_string)
    UniProtXML = etree.fromstring(UniProtXMLstring)
    print 'Number of entries returned from initial UniProt search:', len(UniProtXML)

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

        UniProt_query_string_domains = UniProtXML.xpath('entry/feature[@type="domain"][match_regex(@description, "%s")]' % query_string_domain_selection, extensions = { (None, 'match_regex'): MSMSeeder.core.xpath_match_regex_case_insensitive })

        UniProt_unique_domain_names = set([domain.get('description') for domain in UniProt_query_string_domains])
        print 'Unique domain names selected by the domain selector \'%s\' during the initial UniProt search:\n%s' % (query_string_domain_selection, UniProt_unique_domain_names)

    # =========
    # Print subset of domains returned following filtering with the UniProt_domain_regex (case sensitive)
    # =========

    if UniProt_domain_regex != None:
        regex_matched_domains = UniProtXML.xpath('entry/feature[@type="domain"][match_regex(@description, "%s")]' % UniProt_domain_regex, extensions = { (None, 'match_regex'): MSMSeeder.core.xpath_match_regex_case_sensitive })

        regex_matched_domains_unique_names = set([domain.get('description') for domain in regex_matched_domains])
        print 'Unique domain names selected after searching with the case-sensitive regex string \'%s\':\n%s' % (UniProt_domain_regex, regex_matched_domains_unique_names)

    # =========
    # Now go through all the UniProt entries, domains, PDBs and chains, do some filtering, and store IDs for the selected templates
    # =========

    all_UniProt_entries = UniProtXML.findall('entry')

    selected_templateIDs = []

    print 'Extracting information from returned UniProt data...'

    for entry in all_UniProt_entries:
        entry_name = entry.find('name').text
        if UniProt_domain_regex != None:
            selected_domains = entry.xpath('feature[@type="domain"][match_regex(@description, "%s")]' % UniProt_domain_regex, extensions = { (None, 'match_regex'): MSMSeeder.core.xpath_match_regex_case_sensitive })
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
                PDB_chain_span_nodes = PDB.findall('property[@type="chains"]')

                for PDB_chain_span_node in PDB_chain_span_nodes:
                    chain_span_string = PDB_chain_span_node.get('value')
                    chain_spans = MSMSeeder.UniProt.parse_uniprot_pdbref_chains(chain_span_string)

                    for chainID in chain_spans.keys():
                        span = chain_spans[chainID]
                        if (span[0] < domain_span[0]+30) & (span[1] > domain_span[1]-30):
                            span_string = '%d_%d' % (span[0], span[1])
                            templateID = '%(domainID)s_%(PDBID)s_%(chainID)s_%(span_string)s' % vars()
                            selected_templateIDs.append(templateID)

    print '%d templates selected.' % len(selected_templateIDs)

    # =========
    # Search for PDB and SIFTS files; download if necessary
    # =========


    # =========
    # Write template IDs and sequences to file
    # =========

    with open(fasta_ofilepath, 'w') as fasta_ofile:
        for ID in selected_templateIDs:
            fasta_ofile.write('>' + ID + '\n')

    # =========
    # Update project metadata file
    # =========

    now = datetime.datetime.utcnow()
    datestamp = now.strftime(MSMSeeder.core.datestamp_format_string)

    template_selection_metadata = {
    'datestamp': datestamp,
    'template-selection-method': 'UniProt',
    'UniProt-query-string': UniProt_query_string
    }
    if UniProt_domain_regex != None:
        template_selection_metadata['UniProt-domain-regex'] = UniProt_domain_regex

    project_metadata.data['template-selection'] = template_selection_metadata
    project_metadata.write()

    print 'Done.'

