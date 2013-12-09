# Retrieve domain templates from SCOP given search terms.
#
# John D. Chodera <choderaj@mskcc.org> - 15 Jan 2013; 8 Mar 2013
#
# PREREQUISITES
#
# * MODELLER
# http://salilab.org/modeller/
 
#==============================================================================
# IMPORTS
#==============================================================================

import os.path

#==============================================================================
# PARAMETERS
#==============================================================================

use_regex = False # if True, query terms are used as regex terms
#scop_search_terms = ['Protein kinase-like'] # list of one or more terms or SCOP IDs to match in SCOP description file
scop_search_terms = ['Protein kinases, catalytic subunit'] # list of one or more terms or SCOP IDs to match in SCOP description file

templates_directory = 'templates/' # Destination for templates
templates_filename = os.path.join(templates_directory, 'templates.seg') # Templates file for MODELLER input
templates_index_filename = os.path.join(templates_directory, 'templates.txt') # text index of template filenames
templates_structures_directory = os.path.join(templates_directory, 'structures') # storage for PDB structure files retrieved from SCOP

#==============================================================================
# SCOP file sources.
#==============================================================================

scop_directory = os.path.join('..','structures','scop')
scop_version = "1.75B" # SCOP version
scop_hierarchy_filename      = os.path.join(scop_directory, "dir.hie.scop.%(scop_version)s.txt" % vars())
scop_classification_filename = os.path.join(scop_directory, "dir.cla.scop.%(scop_version)s.txt" % vars())
scop_description_filename    = os.path.join(scop_directory, "dir.des.scop.%(scop_version)s.txt" % vars())
scop_comments_filename       = os.path.join(scop_directory, "dir.com.scop.%(scop_version)s.txt" % vars())

scop_entity_url = "http://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-1.75B/" # location for direct download of SCOP entity domain files

#==============================================================================
# Build a list of SCOP root nodes that match search terms in description file.
#==============================================================================

print "Searching SCOP for query: %s" % str(scop_search_terms)
if use_regex: 
    print "Using query terms as regular expressions."
else:
    print "Using query terms as substrings."

root_node = { 'nodeid' : 0, 'description' : 'root', 'entity' : '-' }
scop_nodes = { '0' : root_node }
root_nodes = list()
infile = open(scop_description_filename)
for line in infile:
    # Skip header lines.
    if line[0] == '#': continue

    # Process entry.
    elements = line.split()
    node = dict()
    node['nodeid'] = elements[0]
    node['level']  = elements[1]
    node['path']   = elements[2]
    node['entity'] = elements[3]
    node['description'] = " ".join(elements[4:]) # reconstitute description

    # Extract parts we'll use for match.
    nodeid = node['nodeid']
    description = node['description']

    # Store node.
    scop_nodes[nodeid] = node
    
    # See if any search terms match SCOP node description.
    match = False
    for query in scop_search_terms:
        if use_regex:
            import re
            if re.search(query, description):
                match = True
                break
        else:
            if description.find(query) != -1:
                match = True
                break

    # Add node to the list if match.
    if match:
        root_nodes.append(nodeid)

infile.close()
print "Query matched %d root nodes." % len(root_nodes)
for nodeid in root_nodes:
    print "%8s %s" % (nodeid, scop_nodes[nodeid]['description'])

#==============================================================================
# Build a list of all SCOP nodes that are children of any root nodes in our list.
#==============================================================================

print "Searching for all child nodes in SCOP hierarchy..."

# Read SCOP hierarchy.
scop_hierarchy = dict()
infile = open(scop_hierarchy_filename)
for line in infile:
    # Skip header lines.
    if line[0] == '#': continue

    # Process entry.
    elements = line.split()
    nodeid = elements[0]
    parent  = elements[1]
    children = elements[2].split(',')
    if children[0] == '-': children = []

    # Store data.
    scop_nodes[nodeid]['parent'] = parent
    scop_nodes[nodeid]['children'] = children

infile.close()

# Build list of child nodes.
last_size = 0
child_nodes = set(root_nodes)
while len(child_nodes) > last_size:
    last_size = len(child_nodes)
    for nodeid in list(child_nodes):
        for child_nodeid in scop_nodes[nodeid]['children']:
            child_nodes.add(child_nodeid)

print "Found a total of %d child nodes (including parent nodes)." % len(child_nodes)
#for nodeid in child_nodes:
#    print "%8s %s" % (nodeid, scop_nodes[nodeid]['description'])

#==============================================================================
# Read all PDB and chain IDs.
#==============================================================================

print "Retrieving PDB and chain IDs..."

# Read SCOP hierarchy.
infile = open(scop_classification_filename)
for line in infile:
    # Skip header lines.
    if line[0] == '#': continue

    # Process entry.
    elements = line.split()
    entityid = elements[0]
    pdbid = elements[1]
    segment = elements[2]
    tree = elements[3]
    nodeid = elements[4]
    
    elements = segment.split(':')
    chainid = elements[0]

    # Store data.
    scop_nodes[nodeid]['pdbid'] = pdbid
    scop_nodes[nodeid]['segment'] = segment
    scop_nodes[nodeid]['chainid'] = chainid

infile.close()

#==============================================================================
# Build a list of all protein domain entities that are assigned to any child node in our list.
#==============================================================================

print "Creating list of all entities..."

entity_nodeids = list()
for nodeid in list(child_nodes):
    if scop_nodes[nodeid]['entity'] != '-':
        entity_nodeids.append(nodeid)

print "Found a total of %d entities." % len(entity_nodeids)

#==============================================================================
# Retrieve structures for entities that do not already exist.
#==============================================================================

import urllib
import os.path

if not os.path.exists(templates_directory):
    os.mkdir(templates_directory)

if not os.path.exists(templates_structures_directory):
    os.mkdir(templates_structures_directory)

for nodeid in entity_nodeids:
    entity_name = scop_nodes[nodeid]['entity']

    # Construct local filename from entity name.
    filename = os.path.join(templates_structures_directory, entity_name + '.pdb')

    # Retrieve coordinate file if it does not exist.
    if not os.path.exists(filename):
        url = scop_entity_url + '/' + entity_name[2:4] + '/' + entity_name + '.ent'
        print "Retrieving '%s' from %s" % (entity_name, url)
        urllib.urlretrieve(url, filename)
        
#==============================================================================
# Perform quality control check on PDB files.
#==============================================================================

print "Performing quality control checks on PDB files..."

import modeller

modeller.log.none()
env = modeller.environ()

# directories for input atom files
env.io.atom_files_directory = templates_structures_directory

for nodeid in entity_nodeids:
    entity_name = scop_nodes[nodeid]['entity']

    # Construct local filename from entity name.
    filename = os.path.join(templates_structures_directory, entity_name + '.pdb')
    
    # Read structure.
    model = modeller.model(env, file=filename)
    
    # Determine number of atoms and residues.
    nresidues = len(model.residues)
    natoms = len(model.atoms)
    
    # Reject if too few atoms.
    MINIMUM_ATOMS_PER_RESIDUE = 3
    if (natoms < MINIMUM_ATOMS_PER_RESIDUE * nresidues):
        # Remove file.
        os.remove(filename)
        
        # Remove entry.
        entity_nodeids.remove(nodeid)
        
        # Print warning.
        print "Removing template %s : %d residues, but only %d atoms found." % (entity_name, nresidues, natoms)
    
#==============================================================================
# Construct initial alignment file for MODELLER.
#==============================================================================

contents = ""
for nodeid in entity_nodeids:
    entity_name = scop_nodes[nodeid]['entity']
    pdbid = scop_nodes[nodeid]['pdbid']
    chainid = scop_nodes[nodeid]['chainid']

    contents += ">P1;" + entity_name + "\n"
    contents += "structureX:%s:FIRST:%s:LAST:%s::::\n" % (entity_name,chainid,chainid)    
    contents += "*\n"

# Write templates filename.
print "Writing templates file for MODELLER to '%s'..." % templates_filename
outfile = open(templates_filename, 'w')
outfile.write(contents)
outfile.close()

#==============================================================================
# Write an index of these templates.
#==============================================================================

contents = ""
for nodeid in entity_nodeids:
    entity_name = scop_nodes[nodeid]['entity']

    # Append filename to index.
    contents += "%s\n" % entity_name

# Write file.
outfile = open(templates_index_filename, 'w')
outfile.write(contents)
outfile.close()

