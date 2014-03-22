import os, urllib, urllib2, datetime, tempfile, subprocess, shutil
# import choderalab as clab
from lxml import etree

def print_uniprot_xml_comparison(new_xml, old_xml):
    '''
    new_xml: lxml root
    old_xml: lxml root
    prints: various statistics
    returns: number of lines in diff between old_xml and new_xml
    '''

    # Get some basic statistics from the XML documents

    new_entries = new_xml.findall('entry')
    new_nentries = len( new_entries)
    new_nentries_with_pdb = len( new_xml.findall('entry/dbReference[@type="PDB"]') )
    new_nentries_with_disease = len( new_xml.findall('entry/comment[@type="disease"]') )

    old_entries = old_xml.findall('entry')
    old_nentries = len( old_entries )
    old_nentries_with_pdb = len( old_xml.findall('entry/dbReference[@type="PDB"]') )
    old_nentries_with_disease = len( old_xml.findall('entry/comment[@type="disease"]') )

    print 'Number of entries (old XML, new XML): %d, %d' % ( old_nentries, new_nentries )
    print 'Number of entries with PDBs (of any kind) (old XML, new XML): %d, %d' % ( old_nentries_with_pdb, new_nentries_with_pdb )
    print 'Number of entries with disease annotations (old XML, new XML): %d, %d' % ( old_nentries_with_disease, new_nentries_with_disease )

    # Also print any ACs that would be removed or added

    new_ACs = [ x.findtext('accession') for x in new_entries ]
    old_ACs = [ x.findtext('accession') for x in old_entries ]

    added_ACs = [ new_AC for new_AC in new_ACs if new_AC not in old_ACs ]
    removed_ACs = [ old_AC for old_AC in old_ACs if old_AC not in new_ACs ]

    if len(added_ACs) > 0:
        print 'ACs to be added:'
        print '\n'.join(added_ACs)
    elif len(removed_ACs) > 0:
        print 'ACs to be removed:'
        print '\n'.join(removed_ACs)

    # Now save xml documents to tempfiles and use diff to get the number of different lines

    try:
        tempdirpath = tempfile.mkdtemp()
        temp_new_xml_filepath = os.path.join(tempdirpath, 'new.xml')
        temp_old_xml_filepath = os.path.join(tempdirpath, 'old.xml')
        with open(temp_new_xml_filepath, 'w') as temp_new_xml_file:
            temp_new_xml_file.write(etree.tostring(new_xml, pretty_print=True))
        with open(temp_old_xml_filepath, 'w') as temp_old_xml_file:
            temp_old_xml_file.write(etree.tostring(old_xml, pretty_print=True))

        temp_diff_filepath = os.path.join(tempdirpath, 'diff.xml')
        with open(temp_diff_filepath, 'w') as temp_diff_file:
            subprocess.call(['diff', '--ignore-matching-lines=<uniprot', temp_old_xml_filepath, temp_new_xml_filepath], stdout=temp_diff_file)
        with open(temp_diff_filepath, 'r') as temp_diff_file:
            len_diff = len( temp_diff_file.readlines() )

    finally:
        shutil.rmtree(tempdirpath)

    print 'Number of lines in diff between the two XML documents:', len_diff

    return len_diff

# This dict converts family information listed in UniProt in the similarity comments to codes similar to those used in the kinase.com poster
# Note that a UniProt "family" is equivalent to a Manning et al. "group". Also, there are a few additional families annotated in UniProt.
kinase_family_uniprot_similarity_text = {
'AGC Ser/Thr protein kinase family' : 'AGC',
'CAMK Ser/Thr protein kinase family' : 'CAMK',
'CMGC Ser/Thr protein kinase family' : 'CMGC',
'CK1 Ser/Thr protein kinase family' : 'CK1',
'STE Ser/Thr protein kinase family' : 'STE',
'TKL Ser/Thr protein kinase family' : 'TKL',
'Tyr protein kinase family' : 'TK',
'NEK Ser/Thr protein kinase family' : 'NEK',
'RIO-type Ser/Thr kinase family' : 'RIO-type'
}

def parse_uniprot_pdbref_chains(chains_span_str):
    '''
    Examples of pdbref chains entries to be parsed:
    A=65-119             => {'A':[65,119]}
    A/C/E/G=64-121       => {'A':[64,121], 'B':[64,121], 'C':[64,121], 'D':[64,121]}
    A=458-778, B=764-778 => {'A':[458,778],'B':[764,778]}
    '''
    comma_sep = chains_span_str.split(',')
    chains_span = dict()
    for s in comma_sep:
        span = s.split('=')[1]
        begin = int(span.split('-')[0])
        end = int(span.split('-')[1])
        chainids = s.split('=')[0].strip().split('/')
        for c in chainids:
            chains_span[c] = [begin, end]
    return chains_span

def encode_uniprot_query(UniProt_query):
    return UniProt_query.replace(' ', '+').replace(':', '%3A').replace('(', '%28').replace(')', '%29')

def retrieve_uniprot(search_string, maxreadlength=100000000):
    '''
    Searches the UniProt database given a search string, and retrieves an XML
    file, which is returned as a string.
    maxreadlength is the maximum size in bytes which will be read from the website
    (default 100MB)
    Example search string: 'domain:"Protein kinase" AND reviewed:yes'

    The function also removes the xmlns attribute from <uniprot> tag, as this
    makes xpath searching annoying
    '''

    base_url = 'http://www.uniprot.org/uniprot/?query='
    search_string_encoded = encode_uniprot_query(search_string)
    query_url = base_url + search_string_encoded + '&format=xml'
    response = urllib2.urlopen(query_url)
    page = response.read(maxreadlength)
    page = page.replace('xmlns="http://uniprot.org/uniprot" ', '', 1)

    return page

def update_metadata_uniprot_search(datestamp, uniprot_search_filepath):
    '''
    Update the datestamp and filepath stored in external-data/metadata.xml
    '''
    now = datetime.datetime.now()
    datestamp = now.strftime(clab.DB.datestamp_format_string)
    parser = etree.XMLParser(remove_blank_text=True)
    metadata_root = etree.parse(clab.DB.external_data_metadata_filepath, parser).getroot()
    uniprot_node = metadata_root.find('UniProt')
    if uniprot_node == None:
        uniprot_node = etree.SubElement(metadata_root, 'UniProt')
    uniprot_search_node = uniprot_node.find('uniprot_search')
    if uniprot_search_node == None:
        uniprot_search_node = etree.SubElement(uniprot_node, 'uniprot_search')
    uniprot_search_node.set('filepath', uniprot_search_filepath)
    uniprot_search_node.set('datestamp', datestamp)
    with open(clab.DB.external_data_metadata_filepath, 'w') as metadata_file:
        metadata_file.write(etree.tostring(metadata_root, pretty_print=True))

def get_uniprot_mapping(query_data_type, retrieve_data_type, query_data):
    '''
NOTE: one gene_id may return multiple uniprotACs.
Use a list of data to query uniprot ID mapping service
and retrieve other types of ID.
query_data can be a list of strings
e.g. get_uniprot_mapping('gene_id','uniprotAC',gene_ids)

Mapping codes can be found here:
http://www.uniprot.org/faq/28#id_mapping_examples
e.g. ACC+ID (from)
     AC (to)
     ID (to)
     P_ENTREZGENEID (both)
    '''

    url = 'http://www.uniprot.org/mapping/'

    query_params = {
    'from':query_data_type,
    'to':retrieve_data_type,
    'format':'tab',
    'query':query_data,
    'reviewed':'yes'
    }

    url_query = urllib.urlencode(query_params)
    request = urllib2.Request(url,url_query)
    request.add_header('User-Agent', 'Python contact')
    response = urllib2.urlopen(request)
    page = response.read(200000)

    # page is in two-column format, with one header line. Extract the info we need into a list
    retrieved_data = page.split('\n')[1:-1]
    for l in range(len(retrieved_data)):
        retrieved_data[l] = retrieved_data[l].split('\t')[1]

    #print url_query
    return retrieved_data

def retrieve_uniprotACs(query_data):
    '''
    Gathers UniProt primary accession IDs by querying UniProt with a given set of GeneIDs.
    Pass it a list of geneID strings
    Returns a list of uniprotACs
    '''
    
    url = 'http://www.uniprot.org/uniprot/'
    retrieved_data = []
    
    for q in query_data:
        print q
        query_params = {
        'query':'geneid:%s AND reviewed:yes' % q,
        'format':'tab',
        'columns':'id'
        }
        
        url_query = urllib.urlencode(query_params)
        request = urllib2.Request(url,url_query)
        request.add_header('User-Agent', 'Python contact')
        response = urllib2.urlopen(request)
        page = response.read(200000)
        
        # each page has one header line. Strip this and put the remaining info into a list
        page_stripped = page.split('\n')[1:-1]
        if len(page_stripped) > 1:
            raise Exception , 'got more than one uniprotAC back'
        # if all is well, then the first value in page_stripped was the only uniprotAC returned
        retrieved_data.append(page_stripped[0])
        
    return retrieved_data

def query_uniprot_multiple(query_params):
    '''
    Queries UniProt with multiple sets of url params.
    Pass: a list of dicts containing params for urllib.urlencode
    Returns: a list of raw response strings.
    '''
    
    url = 'http://www.uniprot.org/uniprot/'
    retrieved_data = []
    
    for q in query_params:
        
        url_query = urllib.urlencode(q)
        request = urllib2.Request(url,url_query)
        request.add_header('User-Agent', 'Python contact')
        response = urllib2.urlopen(request)
        page = response.read(200000)
        
        # each page has one header line. Strip this and put the remaining info into a list
        #page_stripped = page.split('\n')[1:-1]
        #retrieved_data.append(page_stripped)
        retrieved_data.append(page)
        
    return retrieved_data

def query_uniprot(query_params):
    '''
    Queries UniProt with a set of url params.
    Pass: a dict of params for urllib.urlencode
    Returns: the raw response string
    '''
    
    url = 'http://www.uniprot.org/uniprot/'

    url_query = urllib.urlencode(query_params)
    request = urllib2.Request(url,url_query)
    request.add_header('User-Agent', 'Python contact')
    response = urllib2.urlopen(request)
    page = response.read(200000)

    return page

def retrieve_uniprot_xml(uniprotAC):
    '''
    Retrieves a UniProt entry in .xml format
    Pass it a uniprotAC string
    Returns the .xml file as a string
    '''
    
    url = 'http://www.uniprot.org/uniprot/'
    
    response = urllib2.urlopen(url+uniprotAC+'.xml')
    page = response.read(200000)
    
    xml_file = page
        
    return xml_file

