#!/usr/bin/env python
#
# Gathers protein target data - IDs and sequences.
#
# Daniel L. Parton <partond@mskcc.org> - 11 Mar 2014
#

import argparse
argparser = argparse.ArgumentParser(description='Gather target protein data - IDs and sequences.', formatter_class=argparse.RawTextHelpFormatter)

helpstring_gatherfrom = '''Choose a method for gathering target data.
"TargetExplorerDB" (default): Gather target data from an existing
TargetExplorerDB database, specified via either the project metadata file or
the --DBpath argument.
"UniProt": Gather target data from UniProt with a user-defined query
string.'''
argparser.add_argument('--GatherFrom', type=str, help=helpstring_gatherfrom, choices=['TargetExplorerDB', 'UniProt'], default='TargetExplorerDB')
argparser.add_argument('--DBpath', type=str, help='TargetExplorerDB database path', default=None)
args = argparser.parse_args()

def gather_targets_from_TargetExplorerDB():
    '''Gather protein target data from an existing TargetExplorerDB database.'''

    # =========
    # Parameters
    # =========

    import os, datetime, yaml
    import MSMSeeder
    from lxml import etree
    from collections import OrderedDict

    fasta_ofilepath = os.path.join('targets', 'targets.fa')

    project_metadata_filepath = 'project-data.yaml'

    now = datetime.datetime.utcnow()
    datestamp = now.strftime(MSMSeeder.core.datestamp_format_string)

    # =========
    # Read in project metadata file (will throw an exception if it does not exist)
    # =========

    project_metadata = MSMSeeder.core.ProjectMetadata()
    project_metadata.load(project_metadata_filepath)

    # =========
    # Get TargetExplorerDB database path and parse
    # =========

    # First check for command-line arg
    DB_path = args.DBpath

    # DB_path supplied by command-line takes priority. If not found, check the project metadata file
    if DB_path == None:
        try:
            DB_path = project_metadata.data['target-selection']['TargetExplorerDB-database-path']
        except KeyError:
            pass

    if DB_path == None:
        raise Exception, 'TargetExplorerDB database not found in command-line args or in project metadata file. Cannot continue.'

    # Parse the DB
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

    target_selection_metadata = {
    'datestamp': datestamp,
    'target-selection-method': 'TargetExplorerDB',
    'TargetExplorerDB-database-path': DB_path
    }

    project_metadata.data['target-selection'] = target_selection_metadata
    project_metadata.write()

    print 'Done.'


def gather_targets_from_UniProt():
    '''# Searches UniProt for a set of target proteins with a user-defined
    query string, then saves target IDs and sequences.'''

    raise Exception, 'Not yet implemented.'

    # TODO NOTE use 'ABL1_HUMAN_D0' style for targetIDs in this script. Don't bother with mutants. The gather_targets_from_TargetExplorerDB version of this routine will use whatever targetIDs it finds in the DB, which may include 'ABL1_HUMAN_D0_M0' for example.

    # =========
    # Parameters
    # =========

    # =========
    # Get user-defined UniProt query string
    # =========

    # check for command-line arg

    # and check the project metadata file


if __name__ == '__main__':

    if args.GatherFrom == 'TargetExplorerDB':
        gather_targets_from_TargetExplorerDB()
    if args.GatherFrom == 'UniProt':
        gather_targets_from_UniProt()

