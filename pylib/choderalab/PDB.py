# Methods for retrieving data from the PDB
#
# Daniel L. Parton <partond@mskcc.org> - 16 Feb 2013

#==============================================================================
# IMPORTS
#==============================================================================

import urllib2

#==============================================================================
# METHODS
#==============================================================================

def retrieve_sifts(pdb_id):
    '''Retrieves a SIFTS .xml file, given a PDB ID. Works by modifying the PDBe download URL.
    Also removes annoying namespace stuff.
    '''
    import re
    sifts_download_base_url='http://www.ebi.ac.uk/pdbe-srv/view/files/sifts/'
    url = sifts_download_base_url + pdb_id.lower() + '.xml'
    response = urllib2.urlopen(url)
    sifts_page = response.read(10000000) # Max 10MB
    # Old xmlns info
    sifts_page = sifts_page.replace('\nxmlns="http://www.efamily.org.uk/xml/efamily/2004/08/14/eFamily.xsd"\n xmlns:align="http://www.efamily.org.uk/xml/efamily/2004/08/14/alignment.xsd"\n xmlns:data="http://www.efamily.org.uk/xml/data/2004/08/14/dataTypes.xsd"\n xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\nxmlns:dc="http://purl.org/dc/elements/1.1/"\nxmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\nxsi:schemaLocation="http://www.efamily.org.uk/xml/efamily/2004/08/14/eFamily.xsd http://www.efamily.org.uk/xml/efamily/2004/08/14/eFamily.xsd">', '>')
    #sifts_page = sifts_page.replace('align:','')
    #sifts_page = sifts_page.replace('data:','')
    #sifts_page = sifts_page.replace('rdf:','')
    #sifts_page = sifts_page.replace('dc:','')
    #sifts_page = sifts_page.replace('xsi:','')

    # New xmlns info - doesn't work on all files for some reason
    #sifts_page = sifts_page.replace(' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:align="http://www.ebi.ac.uk/pdbe/docs/sifts/alignment.xsd" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:data="http://www.ebi.ac.uk/pdbe/docs/sifts/dataTypes.xsd" xmlns="http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd"', '')
    #sifts_page = sifts_page.replace('schemaLocation="http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd"', '')

    # Ok, I'm just removing all attribs from the entry tag, and the entire rdf tag and contents...
    sifts_page_processed = ''
    skip_rdf_tag_flag = False
    for line in sifts_page.splitlines():
        if line[0:6] == '<entry':
            sifts_page_processed += '<entry>' + '\n'
        elif line[0:7] == '  <rdf:':
            skip_rdf_tag_flag = True
            pass
        elif line[0:8] == '  </rdf:':
            skip_rdf_tag_flag = False
            pass
        else:
            if skip_rdf_tag_flag:
                continue
            sifts_page_processed += line + '\n'
    return sifts_page_processed

def retrieve_pdb(pdb_id,compressed='no'):
    '''Retrieves a PDB file, given a PDB ID. Works by modifying the PDB download URL.
    '''
    pdb_download_base_url='http://www.rcsb.org/pdb/files/'
    url = pdb_download_base_url + pdb_id + '.pdb'
    if compressed == 'yes':
        url += '.gz'
    response = urllib2.urlopen(url)
    pdb_file = response.read(10000000) # Max 10MB
    return pdb_file

