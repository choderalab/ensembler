#!/usr/bin/env python
#
# Searches a user-defined TargetExplorerDB database for a set of target proteins, then saves target IDs and sequences.
#
# Daniel L. Parton <partond@mskcc.org> - 11 Mar 2014
#

# =========
# Parameters
# =========

import os, datetime, json
import MSMSeeder
from lxml import etree
from collections import OrderedDict

fasta_ofilepath = os.path.join('targets', 'targets.fa')

project_metadata_filepath = 'project-data.json'

now = datetime.datetime.utcnow()
datestamp = now.strftime(MSMSeeder.core.datestamp_format_string)

# =========
# Get TargetExplorerDB database path
# =========

# check for command-line arg
# use argparse

# and check the project metadata file

# if query strings are found in both, request user to decide which to use

# XXX temporary
DB_path = '/Users/partond/dev/TargetExplorerDB/database/database.xml'
DB_root = etree.parse(DB_path).getroot()


# update project metadata file
# 'target-selection'
# # 'target-selection-method' node indicating that targets have been selected from a TargetExplorerDB database (rather than from a direct search of UniProt)
# # 'TargetExplorerDB-database-path' - the path to the database
# # 'datestamp'

try:
    with open(project_metadata_filepath, 'r') as project_metadata_file:
        project_metadata = json.load(project_metadata_file)
except IOError as e:
    if e.errno == 2:
        print 'ERROR: Project metadata file "%s" not found. Perhaps you are not in the project top-level directory?' % project_metadata_filepath
    raise

if type(project_metadata) != dict:
    raise Exception, 'Something wrong with the project metadata file. Contained data was not loaded as a dict.'

target_selection_metadata = {
'datestamp': datestamp,
'target-selection-method': 'TargetExplorerDB',
'TargetExplorerDB-database-path': DB_path
}

print target_selection_metadata

project_metadata['target-selection'] = target_selection_metadata

with open(project_metadata_filepath, 'w') as project_metadata_file:
    OD_project_metadata = OrderedDict( sorted(project_metadata.items(), key=lambda x: MSMSeeder.core.project_metadata_document_order.index(x[0])) )
    json.dump(OD_project_metadata, project_metadata_file, indent=4)

# # =========
# # Extract target data from database
# # =========
#
# target_domains = DB_root.findall('entry/UniProt/domains/domain[@targetID]')
#
# # =========
# # Write target data to FASTA file
# # =========
#
# with open(fasta_ofilepath, 'w') as fasta_ofile:
#     for target_domain in target_domains:
#         targetID = target_domain.get('targetID')
#         targetseqnode = target_domain.find('sequence')
#         targetseq = targetseqnode.text.strip()
#
#         target_fasta_string = '>%s\n%s\n' % (targetID, targetseq)
#
#         fasta_ofile.write(target_fasta_string)
#
