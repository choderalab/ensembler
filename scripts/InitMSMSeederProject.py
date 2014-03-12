#!/usr/bin/env python
#
# Initialize project directory by creating necessary subdirectories and a project metadata file.
#
# Daniel L. Parton <partond@mskcc.org> - 11 Mar 2014

import os

# =========
# parameters
# =========

project_dirnames = ['targets', 'templates', 'models', 'packaged-models']
project_metadata_filepath = 'project-data.yaml'

# =========
# initialize project
# =========

for dirname in project_dirnames:
    try:
        os.mkdir(dirname)
        print 'Created directory "%s"' % dirname
    except OSError as e:
        # If directory already exists, e.errno will be set to int 17, and e.strerror will be set to int 'File exists'
        if e.errno == 17:
            print 'Directory "%s" already exists - will not overwrite' % dirname
        else:
            raise
    

if not os.path.exists(project_metadata_filepath):
    with open(project_metadata_filepath, 'w') as project_metadata_file:
        project_metadata_file.write('--- #Project metadata\n')
    print 'Created project metadata file "%s"' % project_metadata_filepath
else:
    print 'Project metadata file "%s" already exists - will not overwrite' % project_metadata_filepath

