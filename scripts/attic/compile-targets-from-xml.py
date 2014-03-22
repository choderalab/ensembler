# Compile sequence data for targets for comparative modeling.
#
# John D. Chodera <choderaj@mskcc.org> - 15 Jan 2013; 8 Mar 2013; 13 Apr 2013 (DLP)
#
# Uses Danny Parton's XML database file constructed from Uniprot as a source
# for sequences.
#
# XXX IMPORTANT: overriding Abl1 sequence so that it includes all of helix I
# (up to residue 513). Eventually will come up with a better automated method
# for determining domain boundaries
#
# PREREQUISITES
#
# * lxml XML parsing library http://lxml.de/
#
# TODO
# 
# * Create a separate directory for each target sequence (with separate
# exceptions.out and sequence identity files) * Split out template sequence
# extraction into a separate script.

import sys, os
from lxml import etree
from choderalab.core import sequnwrap, seqwrap

# PARAMETERS

# Reject sequences below this length.
min_sequence_length = 10

# Print verbose debugging information during execution.
verbose = False 

#
# GET ABSOLUTE PATHS
#

# Input database.
database_filename = '../database/database/kinDB-complete.xml'

# Input directories.
sequence_directory = os.path.abspath("../sequences")

# Output directories.
targets_directory = os.path.abspath("targets")

# Output files.
targets_index_filename = os.path.join(targets_directory, 'targets.txt') # list of target names
targets_sequences_filename = os.path.join(targets_directory, 'targets.fa') # target sequences in FASTA format

# 
# CREATE TARGETS DIRECTORY
#

if not os.path.exists(targets_directory):
    os.mkdir(targets_directory)

#
# BUILD TARGET SEQUENCE LIST IN FASTA FORMAT
#

contents = '' # target sequence file

# Parse targets in sequence database.
parser = etree.XMLParser(remove_blank_text=True)
kinDB = etree.parse(database_filename, parser).getroot()
nkinases = len(kinDB)

target_ids = list()
for k in range(nkinases):
    try:
        # Extract uniprot entry.
        kuniprot = kinDB[k].find('uniprot')
    
        # Extract identity information.
        AC = kuniprot.get('AC')
        entry_name = kuniprot.get('entry_name')

        # Extract PK domains
        pk_domains = kuniprot.findall('pk_domain')
        for pk_domain in pk_domains:
            target_id = pk_domain.get('kinDB_id')
            pk_domain_sequence = sequnwrap(pk_domain.findtext('sequence'))
            len_pk_domain = int(pk_domain.get('length'))

            # XXX XXX XXX IMPORTANT: overriding Abl1 sequence so that it includes all of helix I (up to residue 513). Eventually will come up with a better automated method for determining domain boundaries
            if target_id == 'ABL1_HUMAN_P00519_PK0':
                pk_domain_sequence = sequnwrap(pk_domain.getparent().findtext('sequence'))[241:513]
                len_pk_domain = len(pk_domain_sequence)
        
            # Set target name.
            target_ids.append(target_id)

            # Write alignment file entry.
            contents +=  ">%s\n" % target_id
            contents += seqwrap(pk_domain_sequence)

            if verbose:
                print "%24s : %s" % (target_id, pk_domain_sequence)

        # Mutants
        mutants = kinDB[k].findall('mutants/mutant')
        for mutant in mutants:
            # XXX Skipping these for now - don't have a stable system for assigning IDs yet
            continue
            mut_pk_domain_id = mutant.get('pk_domain_id')
            mutated_full_sequence = list( sequnwrap( kuniprot.find('sequence').text ) )
            pk_domain_begin = int( kuniprot.find('pk_domain[@id="%s"]' % mut_pk_domain_id).get('begin') )
            pk_domain_end = int( kuniprot.find('pk_domain[@id="%s"]' % mut_pk_domain_id).get('end') )
            # XXX IMPORTANT: override Abl1 sequence
            if target_id == 'ABL1_HUMAN_P00519_PK0':
                pk_domain_end = 513

            # Each mutant has >= 1 mutant elements - these will be used to modify mutated_full_sequence
            mutations = mutant.findall('mutation')
            for mutation in mutations:
                orig_resname = mutation.text[0]
                mut_resname = mutation.text[-1]
                mut_resid = int(mutation.text[1:-1])
                if mutated_full_sequence[mut_resid - 1] != orig_resname:
                    raise Exception, 'Wild-type residue in mutation entry does not match that in the stored sequence.'
                mutated_full_sequence[mut_resid - 1] = mut_resname

            # Extract the pk_domain region from the mutated sequence
            mutated_pk_domain_sequence = ''.join( mutated_full_sequence[pk_domain_begin - 1 : pk_domain_end] )

            # Set target name.
            target_id = mutant.get('target_id')
            target_ids.append(target_id)

            # Write alignment file entry.
            contents +=  ">%s\n" % target_id
            contents += seqwrap(mutated_pk_domain_sequence)

            if verbose:
                print "%24s : %s" % (target_id, mutated_pk_domain_sequence)

    except Exception as e:
        import traceback
        print traceback.format_exc()
        print str(e)

outfile = open(targets_sequences_filename, 'w')
outfile.write(contents)
outfile.close()

#
# BUILD TARGETS INDEX FILE
#

contents = ""
for target_id in target_ids:
    contents += "%s\n" % target_id
    
outfile = open(targets_index_filename, 'w')
outfile.write(contents)
outfile.close()

if verbose: print "\n%d target sequences written." % len(target_ids)


