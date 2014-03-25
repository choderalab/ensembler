# Write statistics on progress of kinase simulation set-up procedures.
#
# Run from $CHODERALAB/kinome/models
#
# Daniel L. Parton <partond@mskcc.org> - 29 Oct 2013
#

# ==============
# Imports
# ==============

import sys, os, re, datetime
from lxml import etree
from ast import literal_eval
import choderalab as clab

# ==============
# Parameters
# ==============

parser = etree.XMLParser(remove_blank_text=True)
targets = clab.core.parse_target_args(sys.argv, revert_to_all_targets=True)

models_dir = 'models'
projects_dir = 'projects'
targets_filepath = 'targets/targets.txt'
templates_filepath = 'templates/templates.txt'

statsfilepath = 'stats.xml'

if targets == 'all_targets':
    with open(targets_filepath, 'r') as targets_file:
        targets = targets_file.read().splitlines()

# ==============
# Main
# ==============

# Check for existing stats file
if os.path.exists(statsfilepath):
    stats_root = etree.parse(statsfilepath, parser).getroot()
else:
    stats_root = etree.Element('modelling_progress')

ntargets = clab.core.count_lines_in_file(targets_filepath)
ntemplates = clab.core.count_lines_in_file(templates_filepath)

now = datetime.datetime.now()
datestamp = now.strftime('%Y-%m-%d')

stats_root.set('ntargets', str(ntargets))
stats_root.set('ntemplates', str(ntemplates))
stats_root.set('datestamp', datestamp)

for target in targets:
    # Get directory contents from filesystem
    target_models_dir = os.path.join(models_dir, target)
    target_projects_dir = os.path.join(projects_dir, target)
    if not os.path.exists(target_models_dir):
        continue

    print 'Working on target %s...' % target

    # Look for target XML node - create if it is not found
    target_node = stats_root.find(target)
    if target_node == None:
        target_node = etree.SubElement(stats_root, target)

    # Directories in models/target/
    dirs_in_target_models_dir = os.walk(target_models_dir).next()[1] # picks up only dirs within the target dir
    model_dirs = [ dirpath for dirpath in dirs_in_target_models_dir if clab.core.match_kinDB_ID(dirpath) ]

    # Directories in projects/target/
    try:
        dirs_in_target_projects_dir = os.walk(target_projects_dir).next()[1] # picks up only dirs within the target dir
    except StopIteration:
        dirs_in_target_projects_dir = []
    project_model_dirs = [ dirpath for dirpath in dirs_in_target_projects_dir if clab.core.match_kinDB_ID(dirpath) ]

    # Contents in models/target/model
    model_dirs_contents = [ os.walk(os.path.join(target_models_dir, model_dir)).next()[2] for model_dir in model_dirs ]

    # Contents in models/target/model
    project_model_dirs_contents = [ os.walk(os.path.join(target_projects_dir, project_model_dir)).next()[2] for project_model_dir in project_model_dirs ]

    # Count successful modelling runs
    nmodels = len( [ 1 for model_dir_contents in model_dirs_contents if 'model.pdb' in model_dir_contents ] )
    target_node.set('nmodels', str(nmodels))

    # Count uniques
    nuniques = len( [ 1 for model_dir_contents in model_dirs_contents if 'unique_by_clustering' in model_dir_contents ] )
    target_node.set('nuniques', str(nuniques))

    # Count successful implicit simulation runs
    nimplicit = len( [ 1 for model_dir_contents in model_dirs_contents if 'implicit-refined.pdb' in model_dir_contents ] )
    target_node.set('nimplicit', str(nimplicit))

    # Count successful solvation runs
    nsolvated = len( [ 1 for model_dir_contents in model_dirs_contents if 'nwaters.txt' in model_dir_contents ] )
    target_node.set('nsolvated', str(nsolvated))

    # Count successful explicit simulation runs
    nexplicit = len( [ 1 for model_dir_contents in model_dirs_contents if 'explicit-refined.pdb.txt' in model_dir_contents ] )
    target_node.set('nexplicit', str(nexplicit))

    # Count successful packaging runs
    npackaged = len( [ 1 for model_dir_contents in model_dirs_contents if 'state9.xml' in model_dir_contents ] )
    target_node.set('npackaged', str(npackaged))

    # Set datesteamp
    target_node.set('datestamp', datestamp)


# ==============
# Write output
# ==============

with open(statsfilepath, 'w') as statsfile:
    statsfile.write(etree.tostring(stats_root, pretty_print=True))

