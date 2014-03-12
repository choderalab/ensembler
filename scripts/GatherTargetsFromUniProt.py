#!/usr/bin/env python
#
# Searches UniProt for a set of target proteins with a user-defined query string, then saves target IDs and sequences.
#
# Daniel L. Parton <partond@mskcc.org> - 11 Mar 2014
#
# TODO use 'ABL1_HUMAN_D0' style for targetIDs in this script. Don't bother with mutants. GatherTargetsFromTargetExplorerDB.py will use whatever targetIDs it finds in the DB, which may include 'ABL1_HUMAN_D0_M0' for example.

# =========
# Parameters
# =========

import sys

# =========
# Get user-defined UniProt query string
# =========

# check for command-line arg

# and check the project metadata file

# if query strings are found in both, request user to decide which to use

