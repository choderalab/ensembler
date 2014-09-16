# Extract PDB codes for PK domains in UniProt XML document, download PDB and
# SIFTS XML files (if necessary), write out template.fa and template.txt, and
# save template structures.
#
# New UniProt XML document downloaded only if existing one is > 7 days old, or
# if --forcedl flag is used.
#
# XXX IMPORTANT: overriding Abl1 sequence so that it includes all of helix I
# (up to residue 513). Eventually will come up with a better automated method
# for determining domain boundaries
#
# Daniel L. Parton <parton@cbio.mskcc.org> - 7 May 2013
#
# Dependencies: diff (GNU)
#

#==============================================================================
# IMPORTS
#==============================================================================

import sys, datetime, os, gzip
from choderalab.uniprot import retrieve_uniprot, parse_uniprot_pdbref_chains
from choderalab.pdb import retrieve_sifts, retrieve_pdb
from choderalab.core import seqwrap
from lxml import etree
from subprocess import call
from Bio.PDB import to_one_letter_code

#==============================================================================
# PARAMETERS
#==============================================================================

debug = False
overwrite_templates = True

templates_dir = 'templates'
template_structures_dir = os.path.join(templates_dir, 'structures')
structures_dir = os.path.join('..', 'structures')
pdb_dir = os.path.join(structures_dir, 'pdb')
sifts_dir = os.path.join(structures_dir, 'sifts')
database_dir = '../database'
templates_fa_filename = os.path.join(templates_dir, 'templates.fa') # Templates file in FASTA format
templates_index_filename = os.path.join(templates_dir, 'templates.txt') # text index of template filenames
templates_resnums_filename = os.path.join(templates_dir, 'templates-resnums.txt') # Resnums (1-based UniProt numbering) of template residues, in CSV format

uniprot_xml_path = os.path.join(templates_dir, 'uniprot-kinases_with_pdbs.xml')
uniprot_search_string = 'domain:"protein kinase" AND reviewed:yes AND database:pdb'
uniprot_search_string_query = '?query=domain%3A%22protein+kinase%22+AND+reviewed%3Ayes+AND+database%3Apdb&sort=score&format=xml' # a bit crude, but urllib.urlencode might not be much better
# 1GQ5 is referenced by kinase P16234. The kinase is not in the actual structure. 3O50 and 3O51 have different UniProt ACs in the UniProt entries and SIFTS entries - the SIFTS one is for an unreviewed TrEMBL entry which may be an isoform of the reviewed UniProt entry. Unclear, so ignoring. Plenty of other strucutres for that kinase. 3OG7 is AKAP9-BRAF fusion protein. Fusion proteins are ignored for the moment.
# 1LEW has misannotated UniProt crossrefs in the sifts xml file
ignore_uniprot_pdbs = ['1GQ5', '3O50', '3O51', '3OG7', '1LEW']

if '--forcedl' in sys.argv:
    force_uniprot_download = True
else:
    force_uniprot_download = False

sifts_mappings_path = os.path.join(sifts_dir, 'pdb_chain_uniprot.csv')

ignore_misannotated_sifts = {
'1Q4K':'D'
}  # dict structure is pdbid:chainid

#==============================================================================
# MAIN
#==============================================================================

# =================================
# get sifts structure-level mapping info
# =================================
sifts_mappings_file = open(sifts_mappings_path, 'r')
sifts_mappings_file.readline() # skip first line
sifts_mappings = [ x.strip().split(',') for x in sifts_mappings_file.readlines() ]

if not os.path.exists(templates_dir):
    os.mkdir(templates_dir)

# =================================
# Download UniProt XML document unless it already exists
# =================================
if os.path.exists(uniprot_xml_path):
    print 'UniProt XML document found at:', uniprot_xml_path
# if not, search the UniProt database and save an XML document
else:
    print 'UniProt XML document not found.'
    print 'Retrieving new XML document from UniProt website.'
    new_xml = retrieve_uniprot(uniprot_search_string_query)
    print 'Saving new XML document as:', uniprot_xml_path
    with open(uniprot_xml_path, 'w') as uniprot_xml_file:
        uniprot_xml_file.write(new_xml + '\n')

# =================================
# Parse UniProt kinases XML document
# =================================
print 'Reading UniProt XML document:', uniprot_xml_path
parser = etree.XMLParser(remove_blank_text=True)
uniprot_xml = etree.parse(uniprot_xml_path, parser).getroot()
uniprot_kinases = uniprot_xml.findall('./entry')

# =================================
# Check when the UniProt XML document was retrieved and call retrieve_uniprot() if more than a week old
# =================================
retrieved = uniprot_xml.find('.').attrib['retrieve_date']
retrieved = datetime.datetime.strptime(retrieved, '%Y-%m-%d') # turn into datetime object
now = datetime.datetime.now()
time_elapsed = now - retrieved
print 'UniProt XML document was retrieved: %s (%s days ago)' % (retrieved.strftime('%Y-%m-%d'), time_elapsed.days)
if (time_elapsed.days > 7) or (force_uniprot_download == True):
    if time_elapsed.days > 7:
        print 'UniProt XML document more than 7 days old.'
    elif force_uniprot_download == True:
        print 'Forcing download of new UniProt XML document.'
    print 'Retrieving new XML document from UniProt website.'
    new_xml = retrieve_uniprot(uniprot_search_string_query)
    # First save as a tmp file and carry out a diff with the current XML document
    # There is a python library libdiff, but it is extremely slow compared to the GNU tool
    print 'Conducting diff against current XML document.'
    new_uniprot_xml_path = 'tmp-new-uniprot.xml'
    with open(new_uniprot_xml_path, 'w') as new_uniprot_xml_file:
        new_uniprot_xml_file.write(new_xml + '\n')
    with open('tmp-diff', 'w') as diff_file:
        call(['diff', '--ignore-matching-lines=<uniprot retrieve_date', uniprot_xml_path, new_uniprot_xml_path], stdout=diff_file)
    with open('tmp-diff', 'r') as diff_file:
        diff = diff_file.readlines()
    if len(diff) > 0:
        print 'Differences found:\n==========\n'
        print ''.join(diff)
        print '\n==========\n'
    else:
        print '\n==========\nNo differences found. Continuing anyway.\n==========\n'
    print 'Saving new XML document as:', uniprot_xml_path
    os.rename('tmp-new-uniprot.xml', uniprot_xml_path)
    os.remove('tmp-diff')
    uniprot_xml = etree.parse(uniprot_xml_path, parser).getroot()

nuniprot_kinases = len(uniprot_kinases)
print 'Number of entries in UniProt XML document:', nuniprot_kinases

# =================================
# Iterate through the kinases in the downloaded UniProt XML document, and store relevant info (identifiers, pk_domains, pk_pdbs) in XML etree "templates" (note that this XML etree is only used internally, and is not output to file)
# =================================
templates = etree.Element('templates')
template_ids = []
npdb_chains = 0
pdbids = [] # Just used for checking how many PDB and SIFTS XML files will be required
for k in range(nuniprot_kinases):
    AC = uniprot_kinases[k].findtext('accession')
    entry_name = uniprot_kinases[k].findtext('name')
    sequence = ''.join( uniprot_kinases[k].findtext('sequence[@length]').split() )

    # = UniProt "Protein kinase" domain annotations =
    PK_domains = uniprot_kinases[k].xpath('./feature[@type="domain"][contains(@description,"Protein kinase")]')
    # XXX: SPECIAL CASES
    # These are the entries for which "Protein kinase" domains are known to be not found (case sensitive):
    # MHCKA_DICDI, MYLK_MELGA, TRPM7_MOUSE, TRPM7_RAT
    # These are all alpha-kinases and will be ignored
    if len(PK_domains) < 1:
        print 'Kinase has no "*Protein kinase*" domain annotated. Ignoring.', AC, entry_name
    # Kinases with two PK domains
        continue
    if len(PK_domains) > 1:
        if entry_name in ['E2AK4_HUMAN', 'E2AK4_MOUSE']:
            print 'Kinase contains 2 kinase domains. First is a pseudokinase and will be ignored.', AC, entry_name
            PK_domains.pop(0)
        elif entry_name == 'GCN2_YEAST':
            print 'Kinase contains 2 kinase domain annotations in UniProt, although the presence of two catalytic domains does not seem to be discussed much in literature. Second kinase domain appears to be the active one, although it requires a certain phosphorylation or mutations to be active. PDB structures exist for this domain, but it looks weird, so ignoring.', AC, entry_name
            continue
        elif entry_name in ['JAK1_HUMAN','JAK2_HUMAN','JAK3_HUMAN','JAK2_MOUSE']:
            print 'Kinase contains 2 kinase domains. First is a pseudokinase and will be ignored.', AC, entry_name
            PK_domains.pop(0)
        elif entry_name in ['KS6A1_HUMAN','KS6A2_HUMAN','KS6A3_HUMAN','KS6A4_HUMAN','KS6A5_HUMAN','KS6A6_HUMAN','KS6A3_MOUSE']:
            # Two kinase domains, both active
            pass
        elif entry_name == 'KS6C1_HUMAN':
            print 'Kinase contains 2 kinase domains. First is a pseudokinase and will be ignored.', AC, entry_name
            PK_domains.pop(0)
        elif entry_name == 'OBSCN_HUMAN':
            # Two kinase domains, both active
            pass
        elif entry_name == 'SPEG_HUMAN':
            # Two kinase domains, both active
            pass
        elif entry_name in ['TAF1_HUMAN','TAF1_DROME']:
            # Two kinase domains, both active
            pass
        elif entry_name in ['TYK2_HUMAN','TYK2_MOUSE']:
            # Two kinase domains, both active
            pass
        elif entry_name == 'UNC89_CAEEL':
            print 'Kinase contains 2 kinase domains and is 8081 aa long. Component of C Elegans muscle M-line. Ignoring.', AC, entry_name
            continue
        else:
            raise Exception, 'More than 1 domain found containing "Protein kinase". Please check the following kinase and adjust the script: %s' % entry_name
    # Currently no kinases with PDB structures have a kinase domain annotation which is not equal to "Protein kinase" or "Protein kinase [12]", but just in case...
    if PK_domains[0].attrib['description'] not in ['Protein kinase','Protein kinase 1','Protein kinase 2']:
        raise Exception, 'Unexpected PK domain description. Please adjust the script: %s %s %s' % (AC, entry_name, PK_domains[0].attrib['description'])
    # Now extract the relevant information about the PK domain
    for x_iter,x in enumerate(PK_domains):
        # First calculate the PK domain length and sequence
        pk_description = x.get('description')
        pk_begin = int( x.find('location/begin').attrib['position'] )
        pk_end = int( x.find('location/end').attrib['position'] )
        # XXX XXX XXX IMPORTANT: overriding Abl1 sequence so that it includes all of helix I (up to residue 513). Eventually hoping to come up with a better automated method for determining domain boundaries
        if entry_name == 'ABL1_HUMAN':
            pk_end = 513
        pk_length = pk_end - pk_begin + 1

    # = PDB entries (from UniProt XML) =
    pdbs = uniprot_kinases[k].findall('dbReference[@type="PDB"]')
    npk_pdb_chains = 0
    for p in pdbs:
        # Only keep XRC structures (no NMR or Model)
        if p.find('property[@type="method"]') == None:
            if p.attrib['id'] == '2LV6':
                continue  # 2LV6 has no method listed - it is actually an NMR structure, including only a very short fragment of the kinase, outside the PK domain
        elif p.find('property[@type="method"]').attrib['value'] == 'X-ray':
            pdbid = p.attrib['id']
            if pdbid in ignore_uniprot_pdbs:
                continue
            #resolution = p.find('property[@type="resolution"]').attrib['value']
            chains_span_str = p.find('property[@type="chains"]').attrib['value']
            chains_span = parse_uniprot_pdbref_chains(chains_span_str)
            for c in chains_span.keys():
                if pdbid in ignore_misannotated_sifts.keys():
                    if c == ignore_misannotated_sifts[pdbid]:
                        continue # Some PDB chains are annotated weirdly in the sifts file.
                npdb_chains += 1
                pdb_begin = chains_span[c][0]
                pdb_end = chains_span[c][1]
                # Use the begin and end info to decide if this pdb chain includes the pk_domain. But we will get other sequence info from sifts XML files, using gather-pdb.py
                # Have to check against each PK domain
                for x_iter,x in enumerate(PK_domains):
                    pk_begin = int( x.find('location/begin').attrib['position'] )
                    pk_end = int( x.find('location/end').attrib['position'] )
                    # XXX XXX XXX IMPORTANT: overriding Abl1 sequence so that it includes all of helix I (up to residue 513). Eventually hoping to come up with a better automated method for determining domain boundaries
                    if entry_name == 'ABL1_HUMAN':
                        pk_end = 513
                    pk_length = pk_end - pk_begin + 1
                    if (pdb_begin < pk_begin+30) & (pdb_end > pk_end-30):
                        if pk_length < 350:   # Ignore pk_domains longer than 350 aa
                            template_ids.append( entry_name + '_' + AC + '_PK' + str(x_iter) + '_' + pdbid + '_' + c )
                            template = etree.SubElement(templates, 'template')
                            template.set('template_id', template_ids[-1])
                            template.set('uniprotAC', AC)
                            template.set('pdbid', pdbid)
                            template.set('pk_chainid_uniprot', c)
                            template.set('pk_begin_uniprot', str(pk_begin))
                            template.set('pk_end_uniprot', str(pk_end))
                            template.set('pk_length_uniprot', str(pk_length))
                            npk_pdb_chains += 1
                            if pdbid not in pdbids:
                                pdbids.append(pdbid)

ntemplates = len(template_ids)

# =================================
# Download .pdb and sifts .xml files
# =================================
print '\n= Require %s .pdb and sifts .xml files (will be downloaded if they don\'t exist)=\n' % len(pdbids)
if not os.path.exists(structures_dir):
    print 'Making directory', structures_dir
    os.mkdir(structures_dir)
for t in template_ids:
    # Get the pdbid by splitting the template_id by '_'
    pdbid = t.split('_')[4]
    pdb_filename = os.path.join(pdb_dir, pdbid + '.pdb')
    pdbgz_filename = os.path.join(pdb_dir, pdbid + '.pdb.gz')
    sifts_filename = os.path.join(sifts_dir, pdbid + '.xml.gz')
    if not os.path.exists(pdb_filename):
        print 'Downloading PDB file for:', pdbid
        # Download compressed file, then uncompress
        pdbgz_page = retrieve_pdb(pdbid, compressed='yes')
        with open(pdbgz_filename, 'w') as pdbgz_file:
            pdbgz_file.write(pdbgz_page)    
        with gzip.open(pdbgz_filename, 'rb') as pdbgz_file_decoded:
            with open(pdb_filename, 'w') as pdb_filename:
                pdb_filename.writelines(pdbgz_file_decoded)
        os.remove(pdbgz_filename)
    if not os.path.exists(sifts_filename):
        print 'Downloading sifts file for:', pdbid
        sifts_page = retrieve_sifts(pdbid)
        with gzip.open(sifts_filename, 'wb') as sifts_file: 
            sifts_file.write(sifts_page)    

if not os.path.exists(template_structures_dir):
    os.mkdir(template_structures_dir)

# =================================
# For each template, use the SIFTS XML document to get the PDB resids for observed residues within the pk_domain.
# Then parse the PDB files and output the relevant atoms to template PDB files (in templates/structures).
# =================================
for t in range(ntemplates):
    template = templates[t]
    template_id = template.get('template_id')
    pdbid = template.get('pdbid')
    # DEBUG
    #if pdbid != '4BCP':
    #if pdbid != '4I5C':
    #    continue

    pk_length_uniprot = template.get('pk_length_uniprot')
    pk_begin_uniprot = template.get('pk_begin_uniprot')
    pk_end_uniprot = template.get('pk_end_uniprot')
    pk_chainid_uniprot = template.get('pk_chainid_uniprot')
    if debug:
        print 'Template:', t, template_id, 'pk_domain:', pk_begin_uniprot, pk_end_uniprot, pk_chainid_uniprot,
    else:
        print 'Template:', t, template_id, 'pk_domain:', pk_begin_uniprot, pk_end_uniprot, pk_chainid_uniprot
    # Open and parse SIFTS XML document
    sifts_filename = os.path.join(sifts_dir, pdbid + '.xml.gz')
    with gzip.open(sifts_filename, 'rb') as sifts_file:
        sifts = etree.parse(sifts_file).getroot()

    # XXX: SPECIAL CASES
    # Some PDBs have modified residue(s) (generally TPO or PTR) which have not been given 'PDB modified' annotations in the SIFTS file
    # known cases: 4BCP, 4BCG, 4I5C, 4IVB, 4IAC
    modified_residues = []
    modified_residues.extend( sifts.findall('entity/segment/listResidue/residue[@dbResName="TPO"]') )
    modified_residues.extend( sifts.findall('entity/segment/listResidue/residue[@dbResName="PTR"]') )
    modified_residues.extend( sifts.findall('entity/segment/listResidue/residue[@dbResName="SEP"]') )
    for mr in modified_residues:
        if mr == None:
            continue
        residue_detail_modified = etree.Element('residueDetail')
        residue_detail_modified.set('dbSource','MSD')
        residue_detail_modified.set('property','Annotation')
        residue_detail_modified.text = 'PDB\n          modified'
        mr.append(residue_detail_modified)

    # Try to find the PDB crossref for the first PK domain residue according to the UniProt annotation (only if it is observed)
    pk_begin_pdb_resi = sifts.xpath('entity/segment/listResidue/residue/crossRefDb[@dbSource="UniProt"][@dbResNum="%s"]/../crossRefDb[@dbSource="PDB"][@dbChainId="%s"][not(../residueDetail[contains(text(),"Not_Observed")])]' % (pk_begin_uniprot, pk_chainid_uniprot))
    # If pk_begin_pdb not found, iterate through the UniProt residues until one with a UniProt crossref (and which is observed) is found in the SIFTS file (this should be < 30 aa away, owing to the earlier criterion for picking pk_domains).
    # Note that ensuring we only consider observed residues also means that we don't have to worry (in the case of kinases at least) about PDB resnums such as "324I" which can be a pain, but tend to be non-observed residues.
    if len(pk_begin_pdb_resi) == 0:
        for i in range(1,31):
            pk_begin_pdb_resi = sifts.xpath('entity/segment/listResidue/residue/crossRefDb[@dbSource="UniProt"][@dbResNum="%s"]/../crossRefDb[@dbSource="PDB"][@dbChainId="%s"][not(../residueDetail[contains(text(),"Not_Observed")])]' % ( ( str( int(pk_begin_uniprot) + i ) ), pk_chainid_uniprot) )
            if len(pk_begin_pdb_resi) != 0:
                break
    if len(pk_begin_pdb_resi) == 0:
        if debug:
            print 'pk_begin (or closest residue within 30) not found! Discarding template.'
        template.set('DELETE_ME','')
        template_ids[t] = 'DELETE_ME'
        continue
    pk_begin_pdb = pk_begin_pdb_resi[0].get('dbResNum')

    # Same if pk_end_pdb not found, but iterate in reverse
    pk_end_pdb_resi = sifts.xpath('entity/segment/listResidue/residue/crossRefDb[@dbSource="UniProt"][@dbResNum="%s"]/../crossRefDb[@dbSource="PDB"][@dbChainId="%s"][not(../residueDetail[contains(text(),"Not_Observed")])]' % (pk_end_uniprot, pk_chainid_uniprot))
    if len(pk_end_pdb_resi) == 0:
        for i in range(1,31):
            pk_end_pdb_resi = sifts.xpath('entity/segment/listResidue/residue/crossRefDb[@dbSource="UniProt"][@dbResNum="%s"]/../crossRefDb[@dbSource="PDB"][@dbChainId="%s"][not(../residueDetail[contains(text(),"Not_Observed")])]' % ( ( str( int(pk_end_uniprot) - i ) ), pk_chainid_uniprot) )
            if len(pk_end_pdb_resi) != 0:
                break
    if len(pk_end_pdb_resi) == 0:
        if debug:
            print 'pk_end (or closest residue within 30) not observed! Discarding template.'
        template.set('DELETE_ME','')
        template_ids[t] = 'DELETE_ME'
        continue
    pk_end_pdb = pk_end_pdb_resi[0].get('dbResNum')

    pk_chainid_pdb = sifts.find('entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbResNum="%s"]' % (pk_begin_pdb)).get('dbChainId')
    template.set('pk_chainid_pdb', pk_chainid_pdb)
    if debug:
        print 'SIFTS PDB crossref:', pk_begin_pdb, pk_end_pdb, pk_chainid_pdb

    # Now extract PDB residues which are within the pk_domain bounds, are observed, have a UniProt crossref, and do not have a "PDB modified" or "Conflict" tag.
    # Will also check that at least 70% of the residues within the pk_domain span are observable, and discard the template if not. This is because in some cases (e.g. 3LZB), chains have been added to the PDB entry even though no residues are observable in the structure.
    pk_domain_pdb_length = len( sifts.xpath('entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][@dbResNum >= "%s"][@dbResNum <= "%s"]' % (pk_chainid_pdb, pk_begin_pdb, pk_end_pdb)) )
    pk_domain_pdb_observed = sifts.xpath('entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][@dbResNum >= "%s"][@dbResNum <= "%s"][not(../residueDetail[contains(text(),"Not_Observed")])][../crossRefDb[@dbSource="UniProt"]][not(../residueDetail[contains(text(),"modified")])][not(../residueDetail[contains(text(),"Conflict")])]' % (pk_chainid_pdb, pk_begin_pdb, pk_end_pdb))
    pk_domain_pdb_length_observed = len(pk_domain_pdb_observed)
    ratio_observed = float(pk_domain_pdb_length_observed) / float(pk_domain_pdb_length)
    if debug:
        print pk_domain_pdb_length_observed
    if ratio_observed < 0.7:
        if debug:
            print 'Only %.2f%% residues are observed! Discarding template.' % ratio_observed
        template.set('DELETE_ME','')
        template_ids[t] = 'DELETE_ME'
        continue

    # Generate string of (single-letter) resnames for later outputting to templates.fa, and add as attrib to the templates etree
    # Also make a list of resnums which will be matched against the PDB file
    # And a list of resnums which will be outputted to templates-resnums.txt
    pk_domain_observed_resnames = ''
    pk_domain_observed_pdb_resnums = []
    pk_domain_observed_uniprot_resnums = ''
    for r in pk_domain_pdb_observed:
        residue_uniprot_crossref = r.find('../crossRefDb[@dbSource="UniProt"]')
        pk_domain_observed_resnames += residue_uniprot_crossref.get('dbResName')
        pk_domain_observed_pdb_resnums.append( int( r.get('dbResNum') ) )
        pk_domain_observed_uniprot_resnums += residue_uniprot_crossref.get('dbResNum') + ','
    template.set( 'pk_domain_observed_resnames', pk_domain_observed_resnames )
    template.set( 'pk_domain_observed_uniprot_resnums', pk_domain_observed_uniprot_resnums[0:-1] )  # Leave out final ','

    ### Now parse the PDB file and save the relevant coordinates to template file
    # First build a list of all coords
    pdb_path = os.path.join(pdb_dir, pdbid + '.pdb')
    template_path = os.path.join(template_structures_dir, template_id + '.pdb')
    if overwrite_templates == False and os.path.exists(template_path):
        continue
    with open(pdb_path, 'r') as pdb_file:
        pdb_text = pdb_file.readlines()
    structure_resi_added = {} # used to count how many residues added to the PDB file
    with open(template_path, 'w') as template_file:
        for line in pdb_text:
            if line[0:4] == 'ATOM':
                resnum = int( line[22:26] )
                resname = line[17:20]
                chainid = line[21]
                if chainid == pk_chainid_pdb:
                    if resnum in pk_domain_observed_pdb_resnums:
                        template_file.write(line)
                        structure_resi_added[resnum] = resname
    n_structure_resi_added = len(structure_resi_added)

    # Important check to make sure that the number of residues in the template structure will be the same as that in the sequence
    if debug:
        print 'Number of residues to be added to sequence file:', pk_domain_pdb_length_observed
        print 'Number of residues added to structure file:', n_structure_resi_added
    if n_structure_resi_added != pk_domain_pdb_length_observed:
        structure_one_letter_seq = ''
        for r in sorted(structure_resi_added.keys()):
            structure_one_letter_seq += to_one_letter_code[ structure_resi_added[r] ]
        print 'sequence to be added to sequence file:\n', pk_domain_observed_resnames
        print 'sequence added to structure file:\n', structure_one_letter_seq
        raise Exception('Number of residues to be added to the template structure does not match the number to be added to the sequence file.')

# =================================
# Output templates.txt
# =================================
# First remove templates with DELETE_ME tags (due to low numbers of observed residues)
template_ids_filtered = [ x for x in template_ids if x != 'DELETE_ME' ]
ntemplates_filtered = len(template_ids_filtered)

with open(templates_index_filename, 'w') as templates_index_file:
    for t in template_ids_filtered:
        templates_index_file.write(t + '\n')

# =================================
# Output templates.fa and templates-resnums.txt
# =================================
templates_filtered = templates.xpath('template[not(@DELETE_ME="")]')

with open(templates_fa_filename, 'w') as templates_fa_file:
    with open(templates_resnums_filename, 'w') as templates_resnums_file:
        for t in range(ntemplates_filtered):
            template = templates_filtered[t]
            template_id = template.get('template_id')
            chainid = template.get('pk_chainid_pdb')
            sequence = template.get('pk_domain_observed_resnames')
            resnums = template.get('pk_domain_observed_uniprot_resnums')

            template_header = '>' + template_id + '\n'
            template_fa_string = template_header + seqwrap(sequence)
            template_resnums_string = template_header + resnums + '\n'

            templates_fa_file.write(template_fa_string)
            templates_resnums_file.write(template_resnums_string)

# =================================
# Some stats
# =================================
template_ACs = [ x.get('uniprotAC') for x in templates_filtered ]
nkinases_with_pk_pdb = len(set(template_ACs ))

print 'Total number of pdb chains:', npdb_chains
print '(Number of templates created before filtering: ' + str(ntemplates) + ')'
print 'Total number of templates created:', ntemplates_filtered
print 'Number of kinases with at least one template:', nkinases_with_pk_pdb


