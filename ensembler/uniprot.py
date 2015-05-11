import sys
if sys.version_info > (3, 0):
    from urllib.request import urlopen
else:
    from urllib2 import urlopen
import ensembler
from lxml import etree


def query_uniprot(search_string, maxreadlength=100000000):
    """Searches the UniProt database given a search string, and retrieves an XML
    file, which is returned as a string.
    maxreadlength is the maximum size in bytes which will be read from the website
    (default 100MB)
    Example search string: 'domain:"Protein kinase" AND reviewed:yes'

    The function also removes the xmlns attribute from <uniprot> tag, as this
    makes xpath searching annoying
    """
    base_url = 'http://www.uniprot.org/uniprot/?query='
    search_string_encoded = ensembler.core.encode_url_query(search_string.replace('=', ':'))
    query_url = base_url + search_string_encoded + '&format=xml'
    response = urlopen(query_url)
    page = response.read(maxreadlength)
    page = remove_uniprot_xmlns(page)
    return page


def build_uniprot_query_string_from_acs(acs):
    ac_query_string = ' OR '.join(['acc:%s' % ac for ac in acs])
    return ac_query_string


def get_uniprot_xml(uniprot_query_string, write_to_filepath=None):
    uniprotxmlstring = query_uniprot(uniprot_query_string)
    if len(uniprotxmlstring) == 0:
        raise Exception('UniProt query returned empty string. Query string may have failed to match'
                        ' any UniProt entries, or may have been malformed.')
    if write_to_filepath:
        with open(write_to_filepath, 'w') as ofile:
            ofile.write(uniprotxmlstring)
    parser = etree.XMLParser(huge_tree=True)
    uniprotxml = etree.fromstring(uniprotxmlstring, parser)
    return uniprotxml


def remove_uniprot_xmlns(uniprot_xml_string):
    return uniprot_xml_string.replace('xmlns="http://uniprot.org/uniprot" ', '', 1)


def parse_uniprot_pdbref_chains(chains_span_str):
    """
    Examples of pdbref chains entries to be parsed:
    A=65-119             => {'A':[65,119]}
    A/C/E/G=64-121       => {'A':[64,121], 'B':[64,121], 'C':[64,121], 'D':[64,121]}
    A=458-778, B=764-778 => {'A':[458,778],'B':[764,778]}
    """
    comma_sep = chains_span_str.split(',')
    chains_span = {}
    for s in comma_sep:
        span = s.split('=')[1]
        begin = int(span.split('-')[0])
        end = int(span.split('-')[1])
        chainids = s.split('=')[0].strip().split('/')
        for c in chainids:
            chains_span[c] = [begin, end]
    return chains_span