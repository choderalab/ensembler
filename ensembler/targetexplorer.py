import sys
if sys.version_info > (3, 0):
    from urllib.request import urlopen
else:
    from urllib2 import urlopen
import json


def query_targetexplorer(dbapi_uri, search_string, return_data=None, maxreadlength=10000000):
    """
    Queries a TargetExplorer DB API database given the URI and a search string,
    and returns data as a JSON string.
    maxreadlength is the maximum size in bytes which will be read from the website
    (default 10MB)
    The search string uses SQLAlchemy syntax and standard TargetExplorer
    frontend data fields.
    Example: 'species="Human"'
    Or to select all domains in the DB: ''
    If full_seqs=True, the DB API will also return the full-length canonical
    isoform sequences.
    return_data: str e.g. 'seqs' or list e.g. ['domain_seqs', 'seqs']
    """
    import ensembler.core

    if return_data is None:
        return_data = ''
    elif type(return_data) == str:
        return_data = [return_data]

    base_uri = dbapi_uri + '/search?query='
    search_string_encoded = ensembler.core.encode_url_query(search_string)
    return_data_string = '&return=' + ','.join(return_data)
    query_uri = base_uri + search_string_encoded + return_data_string
    response = urlopen(query_uri)
    page = response.read(maxreadlength)
    return page


def get_targetexplorer_json(dbapi_uri, search_string, return_data=None):
    """
    :param dbapi_uri:
    :param search_string:
    :param return_data: str e.g. 'seqs' or list e.g. ['domain_seqs', 'seqs']
    :return:
    """
    targetexplorer_jsonstr = query_targetexplorer(dbapi_uri, search_string, return_data=return_data)
    return json.loads(targetexplorer_jsonstr)


def get_targetexplorer_metadata(dbapi_uri, maxreadlength=100000):
    """
    Gets metadata for a TargetExplorer DB, via the network API.
    Metadata is returned as a JSON string.
    maxreadlength is the maximum size in bytes which will be read from the website
    (default 100kB)
    """

    full_uri = dbapi_uri + '/get_metadata'
    response = urlopen(full_uri)
    page = response.read(maxreadlength)

    return json.loads(page)