# Build comparative models of each kinase domain from all kinase domain structures in the RCSB using MODELLER.
# This version uses clustal-omega to generate pairwise sequence alignments, rather than the MODELLER automodel.
#
# John D. Chodera <choderaj@mskcc.org> - 15 Jan 2013; 7 Mar 2013; DLP 22 Jun 2013
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

# PARAMETERS

import sys
from ast import literal_eval
import choderalab

# Process only these targets, if specified.
# e.g. -targets '["SRC_HUMAN_PK0_P12931", "ABL1_HUMAN_PK0_P00519"]'
try:
    process_only_these_targets = literal_eval( sys.argv[ sys.argv.index('-targets') + 1 ] )
except ValueError:
    process_only_these_targets = False

if process_only_these_targets:
    print 'Processing only these targets:'
    print process_only_these_targets

process_only_these_templates = False
#process_only_these_templates = ['SRC_HUMAN_P12931_PK0_2H8H_A', 'HCK_HUMAN_P08631_PK0_1AD5_A', 'ABL1_HUMAN_P00519_PK0_2HZ0_A', 'ABL1_HUMAN_P00519_PK0_2HIW_A', 'CDK2_HUMAN_P24941_PK0_1KE7_A']

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

# templates and targets
targets_filename = os.path.join(targets_directory, 'targets.fa')
templates_filename = os.path.join(templates_directory, 'templates.fa')

explicit_hydrogens = False

#
# MODEL ALL SEQUENCES ONTO ALL STRUCTURES
#

def build_model(args):
    """
    Generate a comparative model with MODELLER.
    
    NOTES
    
    Assumes this process has shared filesystem access.
    
    """
    try:
        (target, target_seq, template, template_seq, directories) = args

        # Output.
        print "-------------------------------------------------------------------------"
        print "Modelling '%s' => '%s'" % (target, template)
        print "-------------------------------------------------------------------------"

        # Unpack directories.
        models_directory = directories['models_directory']
        templates_directory = directories['templates_directory']
        clustal_binary = directories['clustal_binary']

        # Construct output filepath and names.
        import os.path
        model_directory = os.path.join(models_directory, target, template)
        output_pdbfilename = os.path.join(model_directory, 'model.pdb')
        seqid_filename = os.path.join(model_directory, 'sequence-identity.txt')
        alignment_filename = os.path.join(model_directory, 'alignment.pir')
        restraint_filename = os.path.join(model_directory, 'restraints.rsr')
        restraint_filename_compressed = os.path.join(model_directory, 'restraints.rsr.gz')

        # Skip model-building if files already exist.
        files_to_check = [model_directory, output_pdbfilename, seqid_filename, alignment_filename, restraint_filename_compressed]
        files_are_present = [os.path.exists(filename) for filename in files_to_check]
        if all(files_are_present):
            text  = "---------------------------------------------------------------------------------\n"
            text += '%s\n' % template
            text += "Output files already exist; were not overwritten.\n"
            return text

        # Initialize MODELLER.
        import modeller
        import modeller.automodel

        modeller.log.none()
        env = modeller.environ()

        # Set directories for input atom files.
        env.io.atom_files_directory = ['.', models_directory, os.path.join(templates_directory, 'structures')]

        if explicit_hydrogens:
            # Use explicit hydrogen atoms.
            env.io.hydrogen = True
            env.libs.topology.read('${LIB}/top.lib')
            env.libs.parameters.read('${LIB}/par.lib')
        
        # Create temporary directory for modeling.
        import os, tempfile
        temporary_directory = tempfile.mkdtemp()
        os.chdir(temporary_directory)

        # Conduct target-template alignment
        seq_ids = [target, template]
        seqs = [target_seq, template_seq]
        aln_seqs = choderalab.align.run_clustalo(seq_ids, seqs)
        target_sequence = aln_seqs[0]
        template_sequence = aln_seqs[1]

        # Write MODELLER PIR file which includes structure file information.
        aligned_filename = 'aligned.pir'
        contents = "Target-template alignment by clustal omega\n"
        contents += ">P1;%s\n" % target
        contents += "sequence:%s:FIRST:@:LAST :@:::-1.00:-1.00\n" % target
        contents += target_sequence + '*\n'
        contents += ">P1;%s\n" % template
        contents += "structureX:%s:FIRST:@:LAST : :undefined:undefined:-1.00:-1.00\n" % template
        contents += template_sequence + '*\n'
        outfile = open('aligned.pir', 'w')
        outfile.write(contents)
        outfile.close()

        # Create an all-atom model.
        a = modeller.automodel.allhmodel(env,
                                         # file with template codes and target sequence
                                         alnfile  = aligned_filename,
                                         # PDB codes of the template
                                         knowns   = template,
                                         # code of the target
                                         sequence = target)
        a.make()                            # do homology modeling

        # Retrieve model.
        model_pdbfilename = a.outputs[0]['name']
        target_model = modeller.model(env, file=model_pdbfilename)

        # Construct directory to place final models in.
        if not os.path.exists(model_directory):
            os.makedirs(model_directory)

        # Write model.
        target_model.write(file=output_pdbfilename)

        # Write sequence identity.
        outfile = open(seqid_filename, 'w')
        outfile.write('%.1f\n' % target_model.seq_id)
        outfile.close()

        # Copy alignment            
        import shutil
        shutil.move('aligned.pir', alignment_filename)

        # Copy restraints.
        import gzip
        with open('%s.rsr' % target, 'r') as rsrfile:
            with gzip.open(restraint_filename_compressed, 'wb') as rsrgzfile:
                rsrgzfile.write(rsrfile.read())

        # Clean up temporary directory.
        os.chdir(current_directory)
        shutil.rmtree(temporary_directory)    

        if os.path.getsize(output_pdbfilename) < 1:
            raise Exception, 'Output PDB file is empty. Could be a filesystem error.'

        text  = "---------------------------------------------------------------------------------\n"
        text += '%s\n' % template
        text += 'Successfully modeled target %s to template %s.\n' % (target, template)
        text += "Sequence identity was %.1f%%.\n" % (target_model.seq_id)
        return text

    except Exception as e:  
        import traceback
        text  = "---------------------------------------------------------------------------------\n"
        text += '%s\n' % template
        text += 'Error modeling target %s to template %s.\n' % (target, template)
        text += str(e) + '\n'
        text += traceback.format_exc()
        return text

if __name__ == "__main__":

    #
    # READ TEMPLATE AND TARGET INDICES AND SEQUENCES
    #

    targets, target_seqs = choderalab.core.parse_fasta_file(targets_filename)
    templates, template_seqs = choderalab.core.parse_fasta_file(templates_filename)
    
    #
    # BUILD MODELS
    # 
    
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

        target_seq = target_seqs[target_iter]

        directories = { 'models_directory' : models_directory, 'templates_directory' : templates_directory, 'clustal_binary' : clustal_binary }
        args_list = list()
        for template_iter, template in enumerate(templates):
            if process_only_these_templates and (template not in process_only_these_templates): continue
            template_seq = template_seqs[template_iter]
            args = (target, target_seq, template, template_seq, directories)
            args_list.append(args)

        from multiprocessing import Pool
        pool = Pool()
        results = pool.map(build_model, args_list)        

        # Create exceptions file.
        log_filename = os.path.join(target_directory, 'models.log')
        outfile = open(log_filename, 'a')
        for result in results:
            outfile.write(result)
        outfile.close()


