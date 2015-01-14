import os
import urllib
import urllib2
import tempfile
import subprocess
import shutil
from core import logger
import ensembler
import ensembler.core
from lxml import etree

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


def query_uniprot(search_string, maxreadlength=100000000):
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
    search_string_encoded = ensembler.core.encode_url_query(search_string.replace('=', ':'))
    query_url = base_url + search_string_encoded + '&format=xml'
    response = urllib2.urlopen(query_url)
    page = response.read(maxreadlength)
    page = remove_uniprot_xmlns(page)
    return page


def build_uniprot_query_string_from_acs(acs):
    ac_query_string = ' OR '.join(['acc:%s' % ac for ac in acs])
    return ac_query_string


def get_uniprot_xml(uniprot_query_string):
    uniprotxmlstring = query_uniprot(uniprot_query_string)
    parser = etree.XMLParser(huge_tree=True)
    uniprotxml = etree.fromstring(uniprotxmlstring, parser)
    return uniprotxml


def remove_uniprot_xmlns(uniprot_xml_string):
    return uniprot_xml_string.replace('xmlns="http://uniprot.org/uniprot" ', '', 1)


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