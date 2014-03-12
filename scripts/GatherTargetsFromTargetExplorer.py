#!/usr/bin/env python
#
# Searches a user-defined TargetExplorerDB database for a set of target proteins, then saves target IDs and sequences.
#
# Daniel L. Parton <partond@mskcc.org> - 11 Mar 2014
#

# =========
# Parameters
# =========

import os
from lxml import etree

fasta_ofilepath = os.path.join('targets', 'targets.fa')

# =========
# Get TargetExplorerDB database path
# =========

# check for command-line arg
# use argparse

# and check the project metadata file

# if query strings are found in both, request user to decide which to use

# TODO update project metadata file with:
# string indicating that targets have been selected from a TargetExplorerDB database (rather than from a direct search of UniProt)
# path to the database
# datestamp for target retrieval


# XXX temporary
DB_path = '/Users/partond/dev/TargetExplorerDB/database/database.xml'
DB_root = etree.parse(DB_path).getroot()

# =========
# Extract target data from database
# =========

target_domains = DB_root.findall('entry/UniProt/domains/domain[@targetID]')

# =========
# Write target data to FASTA file
# =========

with open(fasta_ofilepath, 'w') as fasta_ofile:
    for target_domain in target_domains:
        targetID = target_domain.get('targetID')
        targetseqnode = target_domain.find('sequence')
        targetseq = targetseqnode.text.strip()

        target_fasta_string = '>%s\n%s\n' % (targetID, targetseq)

        fasta_ofile.write(target_fasta_string)

