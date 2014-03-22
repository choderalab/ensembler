# Subject models to explicit solvent simulation.
#
# John D. Chodera <choderaj@mskcc.org> - 17 Feb 2013
#
# PREREQUISITES
#
# * OpenMM
# http://simtk.org/home/openmm
#
# TODO
# * use ff99sb-ildn-star
# * trim all systems to have the same number of waters?

# PARAMETERS

import sys
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

verbose = True

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
print '%d target sequences' % len(targets)
#print "targets:"
#print targets

templates_index_filename = os.path.join(templates_directory, 'templates.txt')
infile = open(templates_index_filename, 'r')
templates = [ line.strip() for line in infile ]
infile.close()
print '%d template structures' % len(templates)
#print "templates:"
#print templates

#
# DETERMINE NWATERS
#

original_directory = os.getcwd()

nwaters_list = list()

for target in targets:
    
    # Process only specified targets if directed.
    if process_only_these_targets and (target not in process_only_these_targets): continue

    target_directory = os.path.join(models_directory, target)
    if not os.path.exists(target_directory): continue

    if verbose: print "Determining number of waters in each system from target '%s'..." % target

    # Process all templates.
    for template in templates:

        model_directory = os.path.join(models_directory, target, template)
        if not os.path.exists(model_directory): continue

        try:
            nwaters_filename = os.path.join(model_directory, 'nwaters.txt')
            infile = open(nwaters_filename, 'r')
            line = infile.readline()
            nwaters = int(line)
            infile.close()

            nwaters_list.append(nwaters)
        
        except Exception as e:
            pass

    import numpy
    nwaters_array = numpy.array(nwaters_list)
    nwaters_array.sort()
    
    filename = os.path.join(target_directory, 'nwaters-list.txt')
    outfile = open(filename, 'w')
    for nwaters in nwaters_array:
        outfile.write('%12d\n' % nwaters)
    outfile.close()

    # display statistics
    index68 = int((len(nwaters_array) - 1) * 0.68)
    index95 = int((len(nwaters_array) - 1) * 0.95)
    print "min = %d, max = %d, mean = %.1f, 68%% = %.0f, 95%% = %.0f\n" % (nwaters_array.min(), nwaters_array.max(), nwaters_array.mean(), nwaters_array[index68], nwaters_array[index95])

    filename = os.path.join(target_directory, 'nwaters-max.txt')
    outfile = open(filename, 'w')
    outfile.write('%d\n' % nwaters_array.max())
    outfile.close()

    # Use 95th percentile.
    filename = os.path.join(target_directory, 'nwaters-use.txt')
    outfile = open(filename, 'w')
    outfile.write('%d\n' % nwaters_array[index68])
    outfile.close()
    

