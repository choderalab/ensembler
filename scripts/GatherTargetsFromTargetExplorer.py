#!/usr/bin/env python
#
# Searches a user-defined TargetExplorerDB database for a set of target proteins, then saves target IDs and sequences.
#
# Daniel L. Parton <partond@mskcc.org> - 11 Mar 2014
#

# =========
# Parameters
# =========

import os, datetime, yaml, argparse
import MSMSeeder
from lxml import etree
from collections import OrderedDict

fasta_ofilepath = os.path.join('targets', 'targets.fa')

project_metadata_filepath = 'project-data.yaml'

now = datetime.datetime.utcnow()
datestamp = now.strftime(MSMSeeder.core.datestamp_format_string)

argparser = argparse.ArgumentParser(description='Gather targets from an existing TargetExplorerDB database.')
argparser.add_argument('--DBpath', type=str, help='TargetExplorerDB database path', default=None)
args = argparser.parse_args()

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

