# Align models onto the template-derived model with highest sequence identity for visual inspection of model quality.
#
# John D. Chodera <choderaj@mskcc.org> - 15 Feb 2013
# Edited DLP 19 Jun 2013
#
# PREREQUISITES
#
# * MODELLER
# http://salilab.org/modeller/
#
# TODO
# * Extract templates list from a template sequence file or directories in models/ directory.

# PARAMETERS

import sys, os
from ast import literal_eval

# Process only these targets, if specified.
# e.g. -targets '["SRC_HUMAN_PK0_P12931", "ABL1_HUMAN_PK0_P00519"]'
try:
    process_only_these_targets = literal_eval( sys.argv[ sys.argv.index('-targets') + 1 ] )
except ValueError:
    process_only_these_targets = False

if process_only_these_targets:
    print 'Processing only these targets:'
    print process_only_these_targets

min_sequence_length = 10

# Look for stride program in PATH environment variable
for path in os.environ['PATH'].split(':'):
    if os.path.exists(os.path.join(path, 'stride')):
        stride_path = os.path.join(path, 'stride')
# If stride not found
if 'stride_path' not in locals():
    # Quick bodge in case we are running on OS X with VMD 1.9.1 installed
    if os.path.exists('/Applications/VMD 1.9.1.app/Contents/vmd/stride_MACOSXX86'):
        stride_path = '/Applications/VMD 1.9.1.app/Contents/vmd/stride_MACOSXX86'
    # Otherwise fail
    else:
        raise Exception, 'Stride program not found. Exiting.'

print 'Stride program found at:', stride_path

#
# GET ABSOLUTE PATHS
#

import os.path

# Input files.
targets_directory = os.path.abspath("targets") # target sequences for modeling
templates_directory = os.path.abspath("templates") # template structures for use in modeling
models_directory = os.path.abspath("models")

#
# INITIALIZE MODELLER
#

from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = '.'

# Use explicit hydrogen atoms.
env.io.hydrogen = True
env.libs.topology.read('${LIB}/top.lib')
env.libs.parameters.read('${LIB}/par.lib')

#
# GET ABSOLUTE PATHS
#

import os.path

# Input files.
targets_directory = os.path.abspath("targets") # target sequences for modeling
templates_directory = os.path.abspath("templates") # template structures for use in modeling
models_directory = os.path.abspath("models")

#
# READ TEMPLATE AND TARGET INDICES
#

targets_index_filename = os.path.join(targets_directory, 'targets.txt')
infile = open(targets_index_filename, 'r')
targets = [ line.strip() for line in infile ]
infile.close()
#print "targets:"
#print targets

templates_index_filename = os.path.join(templates_directory, 'templates.txt')
infile = open(templates_index_filename, 'r')
templates = [ line.strip() for line in infile ]
infile.close()
#print "templates:"
#print templates

#
# ALIGN MODELS ONTO REFERENCE
#

import deprecated_commands
import os.path
for target in targets:
    
    # Process only specified targets if directed.
    if process_only_these_targets and (target not in process_only_these_targets): continue

    # Make sure target directory exists.
    target_directory = os.path.join(models_directory, target)
    if not os.path.exists(target_directory): continue

    print "-------------------------------------------------------------------------"
    print "Aligning '%s'" % (target)
    print "-------------------------------------------------------------------------"

    # Compile sequence identities.
    model_pdbfiles = list()
    model_seqids = list()
    for template in templates:

        model_directory = os.path.join(models_directory, target, template)
        if not os.path.exists(model_directory): continue

        seqid_filename = os.path.join(model_directory, 'sequence-identity.txt')
        if not os.path.exists(seqid_filename): continue

        infile = open(seqid_filename, 'r')
        seqid = float(infile.readline().strip())
        pdb_filename = os.path.join(model_directory, 'model.pdb')
        model_pdbfiles.append(pdb_filename)
        model_seqids.append(seqid)

    # Determine template-derived model with highest sequence identity; this is the reference model.
    import numpy
    reference_index = numpy.array(model_seqids).argmax()
    reference_seqid = model_seqids[reference_index]
    reference_pdbfile = model_pdbfiles[reference_index]
    reference_model = model(env, file=reference_pdbfile)

    # Determine residues to use in alignment.
    try:
        # Get secondary structure information from STRIDE.
        import deprecated_commands
        cmd = stride_path + ' ' + os.path.join(models_directory, reference_pdbfile)
        output = deprecated_commands.getoutput(cmd)
        # Parse STRIDE output, retaining only specified secondary structure types. 
        atmsel = selection()
        core_filename = os.path.join(target_directory, 'core.txt')
        outfile = open(core_filename, 'w')
        for line in output.split('\n'):
            if line[0:3] == 'LOC':
                elements = line.split()
                sstype = elements[1]
                segment_start = int(elements[3])
                segment_end = int(elements[6])
                if sstype in ['AlphaHelix', 'Strand']:
                    # Add CA atoms in segment to atom selection.
                    print line
                    atmsel = atmsel | (selection(reference_model.residue_range('%d:' % segment_start, '%d:' % segment_end)) & selection(reference_model).only_atom_types('CA'))
                    outfile.write('%d %d\n' % (segment_start, segment_end))
        outfile.close()
        print atmsel
    except Exception as e:
        print "Could not call STRIDE; using all CA atoms."
        atmsel = selection(reference_model).only_atom_types('CA')

    for template in templates:

        model_directory = os.path.join(models_directory, target, template)
        if not os.path.exists(model_directory): continue

        target_pdbfile = os.path.join(model_directory, 'model.pdb')
        if not os.path.exists(target_pdbfile): continue
        target_model = model(env, file=target_pdbfile)
        
        # Generate alignment between model and reference_model.
        aln = alignment(env)
        aln.append_model(reference_model, 'reference')
        aln.append_model(target_model, 'target')
        
        # Superimpose model on reference.
        r = atmsel.superpose(target_model, aln)
    
        # Overwrite model.
        aligned_pdbfile = os.path.join(model_directory, 'aligned.pdb')
        target_model.write(file=aligned_pdbfile)

