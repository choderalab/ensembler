# Methods for retrieving data from the PDB
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 16 Feb 2013

#==============================================================================
# IMPORTS
#==============================================================================

import urllib2

#==============================================================================
# METHODS
#==============================================================================

def extract_residues_by_resnum(output_file, pdb_input_file, desired_resnums, desired_chainID):
    '''
    Parameters
    ----------
    output_file: string or gzip.file_like
    pdb_input_file: string or gzip.file_like

    Returns
    ----------
    int
        The number of lines extracted.
    '''
    import gzip
    if type(pdb_input_file) in [str, unicode]:
        with gzip.open(pdb_input_file, 'r') as pdb_file:
            pdbtext = pdb_file.readlines()
    else:
        pdbtext = pdb_input_file.readlines()

    # list of resnum strings e.g. ['9', '29', '30B'] must be converted as follows to match the PDB format:
    # ['   9 ', '  29 ', '  30B']
    import re
    desired_resnums = [ '%4s ' % r if re.match('[0-9]', r[-1]) else '%5s' % r for r in desired_resnums ]

    if type(output_file) in [str, unicode]:
        ofile = open(output_file, 'w')
    else:
        ofile = output_file
    try:
        resnums_extracted = {}
        model_index = 0
        for line in pdbtext:
            # For PDBs containing multiple MODELs (e.g. NMR structures), extract data only from the first model, ignore others.
            if line[0:6] == 'MODEL ':
                model_index += 1
                if model_index == 2:
                    break
            if line[0:6] in ['ATOM  ', 'HETATM']:
                resnum = line[22:27]
                chainID = line[21]
                if chainID == desired_chainID:
                    if resnum in desired_resnums:
                        ofile.write(line)
                        resnums_extracted[resnum] = 1
    finally:
        if type(output_file) in [str, unicode]:
            ofile.close()
    return len(resnums_extracted)

def retrieve_sifts(pdb_id):
    '''Retrieves a SIFTS .xml file, given a PDB ID. Works by modifying the PDBe download URL.
    Also removes annoying namespace stuff.
    '''
    import re, gzip, StringIO
    sifts_download_base_url='ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/'
    url = sifts_download_base_url + pdb_id.lower() + '.xml.gz'
    response = urllib2.urlopen(url)

    sifts_page = response.read(100000000) # Max 100MB
    # Decompress string
    sifts_page = gzip.GzipFile(fileobj=StringIO.StringIO(sifts_page)).read()

    # Removing all attribs from the entry tag, and the rdf tag and contents
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

