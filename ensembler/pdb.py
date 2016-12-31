import sys
if sys.version_info > (3, 0):
    from urllib.request import urlopen
    from urllib.error import URLError
    from io import StringIO
else:
    from urllib2 import urlopen, URLError
    from StringIO import StringIO
import gzip
import re
import six

def extract_residues_by_resnum(output_file, pdb_input_file, template):
    """
    Parameters
    ----------
    output_file: string or gzip.file_like
    pdb_input_file: string or gzip.file_like
    """
    if isinstance(pdb_input_file, six.string_types):
        with gzip.open(pdb_input_file, 'r') as pdb_file:
            pdbtext = pdb_file.readlines()
    else:
        pdbtext = pdb_input_file.readlines()

    # list of resnum strings e.g. ['9', '29', '30B'] must be converted as follows to match the PDB format:
    # ['   9 ', '  29 ', '  30B']
    desired_resnums = ['%4s ' % r if re.match('[0-9]', r[-1]) else '%5s' % r for r in template.resolved_pdbresnums]

    if isinstance(output_file, six.string_types):
        ofile = open(output_file, 'w')
    else:
        ofile = output_file
    try:
        resnums_extracted = {}
        model_index = 0
        for bytesline in pdbtext:
            line = bytesline.decode('UTF-8')
            # For PDBs containing multiple MODELs (e.g. NMR structures), extract data only from the first model, ignore others.
            if line[0:6] == 'MODEL ':
                model_index += 1
                if model_index == 2:
                    break
            if line[0:6] in ['ATOM  ', 'HETATM']:
                resnum = line[22:27]
                chainid = line[21]
                if chainid == template.chainid:
                    if resnum in desired_resnums:
                        ofile.write(line)
                        resnums_extracted[resnum] = 1
    except Exception as e:
        print('Exception detected while extracting ATOM/HETATM records:')
        print(e)
    finally:
        if isinstance(output_file, six.string_types):
            ofile.close()
    if len(resnums_extracted) != len(desired_resnums):
        raise Exception(
                'Number of residues (%d) extracted from PDB (%s) for template (%s) does not match desired number of residues (%d).' % (
                    len(resnums_extracted), template.pdbid, template.templateid, len(desired_resnums)
                )
            )


def retrieve_sifts(pdb_id):
    """Retrieves a SIFTS .xml file, given a PDB ID. Works by modifying the PDBe download URL.
    Also removes annoying namespace stuff.
    """
    sifts_download_base_url='ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/'
    url = sifts_download_base_url + pdb_id.lower() + '.xml.gz'
    try:
        response = urlopen(url)
    except URLError:
        print('ERROR downloading SIFTS file with PDB ID: %s' % pdb_id)
        raise

    sifts_page = response.read(100000000) # Max 100MB
    # Decompress string
    sifts_page = gzip.GzipFile(fileobj=StringIO(sifts_page)).read()

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
    """Retrieves a PDB file, given a PDB ID. Works by modifying the PDB download URL.
    """
    pdb_download_base_url='http://www.rcsb.org/pdb/files/'
    url = pdb_download_base_url + pdb_id + '.pdb'
    if compressed == 'yes':
        url += '.gz'
    response = urlopen(url)
    pdb_file = response.read(10000000) # Max 10MB
    return pdb_file


def extract_uniprot_acs_from_sifts_xml(siftsxml):
    uniprot_crossrefs = siftsxml.findall('entity/segment/listResidue/residue/crossRefDb[@dbSource="UniProt"]')
    uniprot_acs = list(set([uniprot_crossref.get('dbAccessionId') for uniprot_crossref in uniprot_crossrefs]))
    return uniprot_acs
