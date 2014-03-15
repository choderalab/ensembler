#!/usr/bin/env python
#
# Gathers protein template data - IDs, sequences and structures.
#
# Daniel L. Parton <partond@mskcc.org> - 11 Mar 2014
#

# ========
# Command-line arguments
# ========

import argparse
argparser = argparse.ArgumentParser(description='Gather template protein data - IDs, sequences and structures.', formatter_class=argparse.RawTextHelpFormatter)

helpstring_gatherfrom = '''Choose a method for selecting the templates.
"TargetExplorerDB" (default): Select templates from an existing
TargetExplorerDB database, specified via either the project metadata file or
the --DBpath argument.
"UniProt": Select templates from UniProt with a user-defined query
string.'''
argparser.add_argument('--GatherFrom', type=str, help=helpstring_gatherfrom, choices=['TargetExplorerDB', 'UniProt'], default='TargetExplorerDB')
argparser.add_argument('--DBpath', type=str, help='TargetExplorerDB database path', default=None)
argparser.add_argument('--UniProtQuery', type=str, help='UniProt query string', default=None)
argparser.add_argument('--UniProtDomainRegex', type=str, help='Optional regular expression for selecting domains from within UniProt entries', default=None)
args = argparser.parse_args()


# ========
# Definitions
# ========

def gather_templates_from_TargetExplorerDB():
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
    # Get TargetExplorerDB database path and parse
    # =========

    # First check for command-line arg
    DB_path = args.DBpath

    # DB_path supplied by command-line takes priority. If not found, check the project metadata file
    if DB_path == None:
        try:
            DB_path = project_metadata.data['template-selection']['TargetExplorerDB-database-path']
        except KeyError:
            pass

    if DB_path == None:
        raise Exception, 'TargetExplorerDB database not found in command-line args or in project metadata file. Cannot continue.'

    # Parse the DB
    DB_root = etree.parse(DB_path).getroot()

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


def gather_templates_from_UniProt():
    '''# Searches UniProt for a set of template proteins with a user-defined
    query string, then saves IDs, sequences and structures.'''

    # =========
    # Parameters
    # =========

    import os, datetime
    import MSMSeeder
    import MSMSeeder.UniProt
    from lxml import etree

    fasta_ofilepath = os.path.join('templates', 'templates.fa')

    # =========
    # Read in project metadata file
    # =========

    project_metadata = MSMSeeder.core.ProjectMetadata()
    project_metadata.load(MSMSeeder.core.project_metadata_filename)

    # =========
    # Get user-defined UniProt query string and domain-selection regex
    # =========

    # First check for command-line arg
    UniProt_query_string = args.UniProtQuery
    UniProt_domain_regex = args.UniProtDomainRegex

    # Args supplied from command-line take priority. If not found, check the project metadata file
    if UniProt_query_string == None:
        try:
            UniProt_query_string = project_metadata.data['template-selection']['UniProt-query-string']
        except KeyError:
            pass
    if UniProt_domain_regex == None:
        try:
            UniProt_domain_regex = project_metadata.data['template-selection']['UniProt-domain-regex']
        except KeyError:
            pass

    if UniProt_query_string == None:
        raise Exception, 'UniProt query string not found in command-line args or in project metadata file. Cannot continue.'
    # UniProt_domain_regex is optional

    print UniProt_query_string

    # =========
    # Make request to UniProt web server
    # =========

    UniProtXMLstring = MSMSeeder.UniProt.retrieve_uniprot(UniProt_query_string)
    UniProtXML = etree.fromstring(UniProtXMLstring)
    print len(UniProtXML)

    # =========
    # Parse returned UniProt data
    # =========

    # =========
    # Search for template PDB and SIFTS files; download if necessary
    # =========

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


if __name__ == '__main__':

    if args.GatherFrom == 'TargetExplorerDB':
        gather_templates_from_TargetExplorerDB()

    if args.GatherFrom == 'UniProt':
        gather_templates_from_UniProt()

