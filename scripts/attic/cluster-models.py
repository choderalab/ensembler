# RMSD-based clustering to filter out models with high structural similarity
#
# Daniel L. Parton <partond@mskcc.org> - 22 Feb 2013
# Edits by John D. Chodera <choderaj@mskcc.org> - 17 Mar 2013
#
# PREREQUISITES
#
# * MDAnalysis
# code.google.com/p/mdanalysis

# =============================
# IMPORTS
# =============================

import os,sys,glob,shutil
from ast import literal_eval
import MDAnalysis as MDA
from MDAnalysis.analysis import rms
from MDAnalysis.coordinates import DCD
from pylab import *

# =============================
# PARAMETERS
# =============================

# Process only these targets, if specified.
# e.g. -targets '["SRC_HUMAN_PK0_P12931", "ABL1_HUMAN_PK0_P00519"]'
try:
    process_only_these_targets = literal_eval( sys.argv[ sys.argv.index('-targets') + 1 ] )
except ValueError:
    process_only_these_targets = False

if process_only_these_targets:
    print 'Processing only these targets:'
    print process_only_these_targets

targets_directory = os.path.abspath("targets")
templates_directory = os.path.abspath("templates")
models_directory = os.path.abspath("models")

verbose = True

cutoff = 0.6  # Cutoff for RMSD clustering (Angstroms)

# Option to output min_rmsd.dat and min_rmsd.png - the distribution of nearest neighbor rmsds
# NOTE: plotting output may fail if this script is run remotely without X forwarding
min_rmsd_option = False

# =============================
# READ TEMPLATE AND TARGET INDICES
# =============================

targets_index_filename = os.path.join(targets_directory, 'targets.txt')
infile = open(targets_index_filename, 'r')
targets = [ line.strip() for line in infile ]
infile.close()
print 'number of targets:', len(targets)

templates_index_filename = os.path.join(templates_directory, 'templates.txt')
infile = open(templates_index_filename, 'r')
templates = [ line.strip() for line in infile ]
infile.close()
print 'number of templates:', len(templates)

for target in targets:
    # Process only specified targets if directed.
    if process_only_these_targets and (target not in process_only_these_targets): continue

    model_directory = os.path.join(models_directory, target)
    if not os.path.exists(model_directory): continue

    # =============================
    # Combine the kinase models into a temporary DCD trajectory
    # =============================

    if verbose: print 'Compiling models into a temporary DCD trajectory...'

    # Create a temporary DCD file compiling all models.
    import os, tempfile
    temporary_directory = tempfile.mkdtemp()

    # Create DCD file.

    # Create first MDA Universe just to get the number of atoms
    dcd_filename = os.path.join(models_directory,target,'tmp-models.dcd')
    u = MDA.Universe(os.path.join(models_directory,target,templates[0],'model.pdb'))
    dcdw = DCD.DCDWriter(dcd_filename, u.atoms.numberOfAtoms())
    valid_templates = list()
    for (t, template) in enumerate(templates): 
        print 'Writing %s model %d of %d\r' % (target, t+1 , len(templates)),
        model_pdb_filename = os.path.join(models_directory,target,template,'model.pdb')
        if os.path.exists(model_pdb_filename):
            u.load_new(model_pdb_filename)
            dcdw.write(u)
            # Record list of valid templates.
            valid_templates.append(template)
    dcdw.close()
    if verbose: print "%d valid models (out of %d templates) compiled" % (len(valid_templates), len(templates))

    # =============================
    # CLUSTERING
    # =============================

    print 'Conducting clustering...'

    # Remove any existing unique_by_clustering files
    for f in glob.glob( os.path.join(models_directory,target)+'/*_PK_*/unique_by_clustering' ):
        os.unlink(f)

    # Take first pdb to get topol info
    pdb_path = os.path.join(models_directory,target,valid_templates[0],'model.pdb')
    # Then load the MDAnalysis Universe
    u = MDA.Universe(pdb_path,dcd_filename)

    # Option to output min_rmsd.dat, and to write the template and its
    # nearest neighbor to a DCD.
    if min_rmsd_option:
        min_rmsd_dcd_path = os.path.join(models_directory,target,'min_rmsd_pairs.dcd')
        dcdw = DCD.DCDWriter(min_rmsd_dcd_path, u.atoms.numberOfAtoms())

    # Each template will be added to the list uniques if it is further than
    # 0.2 Angstroms (RMSD) from the nearest template.
    uniques=[]
    min_rmsd = []
    for (t, template) in enumerate(valid_templates):
        # Add the first template to the list of uniques
        if t==0:
            uniques.append(template)
            continue

        # Cluster using CA atoms
        rmsd_obj = rms.RMSD(u, select='name CA')
        rmsd_obj.run(start=0,stop=t,ref_frame=t)
        rmsd = rmsd_obj.rmsd.swapaxes(0,1)[2]
        min_rmsd.append( min(rmsd) )

        if min_rmsd[-1] < cutoff:
            if min_rmsd_option:
                dcdw.write(u.trajectory[t]) # write the current template
                # find the index of the nearest neighbor
                min_rmsd_index = list(rmsd).index(min_rmsd[-1])
                dcdw.write(u.trajectory[min_rmsd_index]) # write the nearest neighbor
                u.trajectory[t] # go back to the current template
            continue
        else:
            uniques.append( template )
            # Create a blank file to say this template was found to be unique
            # by clustering
            open( os.path.join(models_directory,target,template,'unique_by_clustering'), 'w').close()

        u.trajectory.rewind()

    uniques_file = open( os.path.join(models_directory,target,'unique-models.txt'), 'w')
    for u in uniques:
        uniques_file.write(u+'\n')
    uniques_file.close()
    print '%d unique models (from original set of %d) using cutoff of %.2f Angstroms' % (len(uniques), len(valid_templates), cutoff)

    if min_rmsd_option:
        # Close the DCD (will only be open at this stage if min_rmsd_option is True)
        dcdw.close()
        # Save the min_rmsd data to a text file
        min_rmsd_path = os.path.join(models_directory,target,'min_rmsd.dat')
        savetxt(min_rmsd_path, min_rmsd)
        # And plot a histogram
        plot_path = os.path.join(models_directory,target,'min_rmsd.png')
        bin_ranges=arange(0,18,0.1)
        # Following matplotlib cli_commands may fail if this script is being run on a remote server without X forwarding
        try:
            hist(min_rmsd, bins=bin_ranges)
            xlim(0,2)
            xticks(arange(0,2.1,0.2))
            yticks(arange(0,121,10))
            xlabel('RMSD ($\AA$)')
            ylabel('population')
            grid()
            savefig(plot_path)
        except:
            pass

    # Clean up temporary directory.
    shutil.rmtree(temporary_directory)    

print 'Done.'

