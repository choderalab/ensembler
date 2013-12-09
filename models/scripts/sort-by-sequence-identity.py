# Compile sorted list of templates by sequence identity.
#
# John D. Chodera <choderaj@mskcc.org> - 17 Mar 2013

# PARAMETERS

import sys
from ast import literal_eval
verbose = True

# Process only these targets, if specified.
# e.g. -targets '["SRC_HUMAN_PK0_P12931", "ABL1_HUMAN_PK0_P00519"]'
try:
    process_only_these_targets = literal_eval( sys.argv[ sys.argv.index('-targets') + 1 ] )
except ValueError:
    process_only_these_targets = False

if process_only_these_targets:
    print 'Processing only these targets:'
    print process_only_these_targets

if __name__ == '__main__':

    #
    # GET ABSOLUTE PATHS
    #

    import os.path

    # Input files.
    targets_directory = os.path.abspath("targets") # target sequences for modeling
    templates_directory = os.path.abspath("templates") # template structures for use in modeling
    models_directory = os.path.abspath("models")

    # Output files.
    projects_directory = os.path.abspath("projects") # FAH projects directory

    #
    # READ TEMPLATE AND TARGET INDICES
    #

    targets_index_filename = os.path.join(targets_directory, 'targets.txt')
    infile = open(targets_index_filename, 'r')
    targets = [ line.strip() for line in infile ]
    infile.close()
    print '%d target sequences' % len(targets)

    templates_index_filename = os.path.join(templates_directory, 'templates.txt')
    infile = open(templates_index_filename, 'r')
    templates = [ line.strip() for line in infile ]
    infile.close()
    print '%d template structures' % len(templates)

    #
    # COMPILE SORTED LIST BY SEQUENCE IDENTITY
    #

    for target in targets:
        
        # Process only specified targets if directed.
        if process_only_these_targets and (target not in process_only_these_targets): continue

        target_directory = os.path.join(models_directory, target)
        if not os.path.exists(target_directory): continue

        print "-------------------------------------------------------------------------"
        print "Compiling template sequence identities for target %s" % (target)
        print "-------------------------------------------------------------------------"

        #
        # BUILD A LIST OF VALID TEMPLATES
        #

        # Process all templates.
        if verbose: print "Building list of valid templates..."
        import os.path
        valid_templates = list()
        for template in templates:
            model_filename = os.path.join(target_directory, template, 'model.pdb')
            if os.path.exists(model_filename):
                valid_templates.append(template)

        nvalid = len(valid_templates)
        if verbose: print "%d valid models found" % nvalid

        #
        # SORT BY SEQUENCE IDENTITY
        #

        if verbose: print "Sorting templates in order of decreasing sequence identity..."
        import numpy
        sequence_identities = numpy.zeros([nvalid], numpy.float32)
        for (template_index, template) in enumerate(valid_templates):
            filename = os.path.join(target_directory, template, 'sequence-identity.txt')
            infile = open(filename, 'r')
            contents = infile.readline().strip()
            infile.close()
            sequence_identity = float(contents)
            sequence_identities[template_index] = sequence_identity
        sorted_indices = numpy.argsort(-sequence_identities)

        #
        # WRITE TEMPLATES SORTED BY SEQUENCE IDENTITY
        # 

        filename = os.path.join(target_directory, 'sequence-identities.txt')
        print filename
        outfile = open(filename, 'w')
        for index in sorted_indices:
            template = valid_templates[index]
            identity = sequence_identities[index]
            outfile.write('%-40s %6.1f\n' % (template, identity))
        outfile.close()

