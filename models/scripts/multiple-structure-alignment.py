# Build a structure-structure multiple sequence alignment of a subset of SCOP kinase domain structures.
#
# John D. Chodera <choderaj@mskcc.org> - 15 Jan 2013; 5 Mar 2013
#
# Based on:
# A sample script for fully automated comparative modeling (Ben Webb)
#
# USAGE
#
# python scripts/build-models.py 
#
# PREREQUISITES
#
# * MODELLER
# http://salilab.org/modeller/
#
# TODO
# 
# * Take list of target names or indices as command-line argument

# PARAMETERS

# Process only these targets, if specified.
#process_only_these_targets = False
process_only_these_targets = ['SRC']

verbose = True

#
# GET ABSOLUTE PATHS
#

import sys, os.path

# Input files.
targets_directory = os.path.abspath("targets") # target sequences for modeling
templates_directory = os.path.abspath("templates") # template structures for use in modeling

# Output files.
models_directory = os.path.abspath("models") # directory to organize models, by target

#
# INITIALIZE MODELLER
#

from modeller import *

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = './:%(models_directory)s:%(templates_directory)s/structures' % vars()

# Use explicit hydrogen atoms.
env.io.hydrogen = True
env.libs.topology.read('${LIB}/top.lib')
env.libs.parameters.read('${LIB}/par.lib')

#
# READ TEMPLATE AND TARGET INDICES
#

targets_index_filename = os.path.join(targets_directory, 'targets.txt')
infile = open(targets_index_filename, 'r')
targets = [ line.strip() for line in infile ]
infile.close()
print "targets:"
print targets

templates_index_filename = os.path.join(templates_directory, 'templates.txt')
infile = open(templates_index_filename, 'r')
templates = [ line.strip() for line in infile ]
infile.close()
print "templates:"
print templates

#aln = alignment(env)
#print "Reading templates..."
#aln.append('%(template_directory)s/templates.seg' % vars())
#templates = [ alignment.code for alignment in aln ]

#
# READ TEMPLATES
#

# Select subset of templates to use for alignment.
nskip = 10
templates_subset = templates[::nskip]

print "Reading templates..."
aln = alignment(env)
templates_unaligned_filename = os.path.join(templates_directory, 'templates.seg')
#aln.append(templates_unaligned_filename)
aln.append(templates_unaligned_filename, align_codes=templates_subset)

#
# CONSTRUCT STRUCTURE-BASED TEMPLATE ALIGNMENT
# 

dendrogram_filename = os.path.join(templates_directory, 'dendrogram.tree')
iteration = 1
for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    if verbose: print "Alignment iteration %d..." % iteration
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50),
               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
               dendrogram_file=dendrogram_filename,
               #alignment_type='progressive', # If 'progresive', the tree is not
               alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
               feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
               improve_alignment=True, fit=True, write_fit=write_fit,
               write_whole_pdb=whole, output='ALIGNMENT QUALITY')
    iteration += 1

templates_aligned_pap_filename = os.path.join(templates_directory, 'templates.pap')
templates_aligned_pir_filename = os.path.join(templates_directory, 'templates.pir')
aln.write(file=templates_aligned_pap_filename, alignment_format='PAP')
aln.write(file=templates_aligned_pir_filename, alignment_format='PIR')

if verbose: print "Assessing alignment quality..."
tree_filename = os.path.join(templates_directory, 'templates.tree')
aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
           rr_file='$(LIB)/as1.sim.mat', overhang=30,
           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
           gap_gap_score=0, gap_residue_score=0, dendrogram_file=tree_filename,
           alignment_type='progressive', feature_weights=[0]*6,
           improve_alignment=False, fit=False, write_fit=True,
           write_whole_pdb=False, output='QUALITY')

if verbose: print "Done."
