# Build comparative models of each kinase domain from all kinase domain structures in the RCSB using MODELLER.
# This version uses clustal-omega to generate pairwise sequence alignments, rather than the MODELLER automodel.
#
# John D. Chodera <choderaj@mskcc.org> - 15 Jan 2013; 7 Mar 2013
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
# * Auto-detect architecture for clustal-omega


# PARAMETERS

# Process only these targets, if specified.
#process_only_these_targets = False
# process_only_these_targets = ['SRC']
process_only_these_targets = ['SRC_HUMAN_PK0_P12931'] # SRC Uniprot AC

#
# GET ABSOLUTE PATHS
#

import os.path
import platform
import sys

# Input files.
targets_directory = os.path.abspath("targets") # target sequences for modeling
templates_directory = os.path.abspath("templates") # template structures for use in modeling

# Output files.
models_directory = os.path.abspath("models") # directory to organize models, by target

# clustal-omega binary
if platform.system() == 'Darwin' and platform.machine() == 'x86_64':
    clustal_binary = os.path.abspath("../external-tools/clustal-omega/clustal-omega-1.0.3-Mac-ppc-x86_64")
elif platform.system() == 'Linux' and platform.machine() == 'x86_64':
    clustal_binary = os.path.abspath("../external-tools/clustal-omega/clustalo-1.1.0-linux-64")
else:
    print 'Operating system does not support included clustal binaries. Exiting...'
    sys.exit()

#
# INITIALIZE MODELLER
#

from modeller import *
from modeller.automodel import *    # Load the automodel class

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

#
# READ TEMPLATE AND TARGET SEQUENCES
#

targets_filename = os.path.join(targets_directory, 'targets.fa')
with open(targets_filename, 'r') as targets_file:
    target_strings = targets_file.read().split('>')
    target_seqs = [ ''.join(target_string.split('\n')[1:]) for target_string in target_strings ]
templates_filename = os.path.join(templates_directory, 'templates.fa')
with open(templates_filename, 'r') as templates_file:
    template_strings = templates_file.read().split('>')
    template_seqs = [ ''.join(template_string.split('\n')[1:]) for template_string in template_strings ]

#
# MODEL ALL SEQUENCES ONTO ALL STRUCTURES
#

# Create exceptions file.
error_filename = os.path.join(models_directory, 'exceptions.out')
outfile = open(error_filename, 'w')
outfile.close()

# Get current working directory.
import os
current_directory = os.getcwd() 

for target_iter, target in enumerate(targets):
    
    # Process only specified targets if directed.
    if process_only_these_targets and (target not in process_only_these_targets): continue

    print "-------------------------------------------------------------------------"
    print "Modelling '%s'" % (target)
    print "-------------------------------------------------------------------------"

    # Create directory for target.
    target_directory = os.path.join(models_directory, target)
    if not os.path.exists(target_directory):
        os.makedirs(target_directory)

    # Process all templates.
    for template_iter, template in enumerate(templates):
        try:
            print "-------------------------------------------------------------------------"
            print "'%s' => '%s'" % (template, target)
            print "-------------------------------------------------------------------------"
            
            # Create temporary directory for modeling.
            import os, tempfile
            temporary_directory = tempfile.mkdtemp()
            print temporary_directory
            os.chdir(temporary_directory)

            # Create template-target alignment.
            target_fa_string = '>' + target + '\n' + target_seq[target_iter] + '\n'
            template_fa_string = '>' + template + '\n' + template_seq[template_iter] + '\n'
            input_fa_string = target_fa_string + template_fa_string

            from subprocess import Popen, PIPE
            import shlex
            command = '%(clustal_binary)s -i - --dealign --outfmt=vienna --force' % vars()
            p = Popen(shlex.split(command), stdin=PIPE, stdout=PIPE)
            stdout = p.communicate(input=input_fa_string)[0]
            print stdout

            # Read target and template sequences.
            infile = open('aligned.vienna', 'r')
            infile.readline()
            target_sequence = infile.readline().rstrip('\n')
            infile.readline()
            template_sequence = infile.readline().rstrip('\n')
            infile.close()

            # Write MODELLER PIR file which includes structure file information.
            aligned_filename = 'aligned.pir'
            print "Creating PIR file..."
            contents = "Target-template alignment by clustal omega\n"
            contents += ">P1;%s\n" % target
            contents += "sequence:SRC:FIRST:@:LAST :@:::-1.00:-1.00\n"
            contents += target_sequence + '*\n'
            contents += ">P1;%s\n" % template
            contents += "structureX:%s:FIRST:@:LAST : :undefined:undefined:-1.00:-1.00\n" % template
            contents += template_sequence + '*\n'
            outfile = open('aligned.pir', 'w')
            outfile.write(contents)
            outfile.close()

            # Create an all-atom model.
            a = allhmodel(env,
                          # file with template codes and target sequence
                          alnfile  = aligned_filename,
                          # PDB codes of the template
                          knowns   = template,
                          # code of the target
                          sequence = target)
            a.make()                            # do homology modeling

            # Retrieve model.
            model_pdbfilename = a.outputs[0]['name']
            target_model = model(env, file=model_pdbfilename)

            # Construct directory to place final models in.
            model_directory = os.path.join(models_directory, target, template)
            if not os.path.exists(model_directory):
                os.makedirs(model_directory)

            # Write model.
            output_pdbfilename = os.path.join(model_directory, 'model.pdb')
            target_model.write(file=output_pdbfilename)

            # Write sequence identity.
            seqid_filename = os.path.join(model_directory, 'sequence-identity.txt')
            outfile = open(seqid_filename, 'w')
            outfile.write('%.1f\n' % target_model.seq_id)
            outfile.close()

            # Copy alignment            
            import shutil
            alignment_filename = os.path.join(model_directory, 'alignment.pir')
            shutil.move('aligned.pir', alignment_filename)

            # Copy restraints.
            restraint_filename = os.path.join(model_directory, 'restraints.rsr')
            commands.getoutput('cp %s.rsr %s; gzip %s' % (target, restraint_filename, restraint_filename))

            # Clean up temporary directory.
            os.chdir(current_directory)
            shutil.rmtree(temporary_directory)

        except Exception as e:  
            # If an Exception has occurred, write a note to the exceptions file and keep running.
            outfile = open(error_filename, 'a')
            outfile.write('---------------------------------------------------------------------------------\n')
            outfile.write('Error modeling target %s to template %s.\n' % (target, template))
            outfile.write(str(e) + '\n')
            outfile.close()


