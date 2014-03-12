#!/usr/bin/env python
#
# Initialize project directory by creating necessary subdirectories and a project metadata file.
#
# Daniel L. Parton <partond@mskcc.org> - 11 Mar 2014

# =========
# Parameters
# =========

import os, datetime, json
import MSMSeeder

project_dirnames = ['targets', 'templates', 'models', 'packaged-models']
project_metadata_filepath = 'project-data.json'

now = datetime.datetime.utcnow()
datestamp = now.strftime(MSMSeeder.core.datestamp_format_string)

# =========
# Create necessary project directories
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
    
# =========
# Create metadata file and add datestamp
# =========

if not os.path.exists(project_metadata_filepath):
    with open(project_metadata_filepath, 'w') as project_metadata_file:
        init_node = {'init' : {'datestamp' : datestamp}}
        json.dump(init_node, project_metadata_file, indent=4)
    print 'Created project metadata file "%s"' % project_metadata_filepath
else:
    print 'Project metadata file "%s" already exists - will not overwrite' % project_metadata_filepath

