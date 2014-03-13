#!/usr/bin/env python
#
# Initialize project directory by creating necessary subdirectories and a project metadata file.
#
# Daniel L. Parton <partond@mskcc.org> - 11 Mar 2014

def init_msmseeder_project(project_toplevel_dir):
    '''Initialize MSMSeeder project within the current directory. Creates
    necessary subdirectories and a project metadata .yaml file.
    '''

    # =========
    # Parameters
    # =========

    import os, datetime, yaml, argparse
    import MSMSeeder

    project_dirnames = ['targets', 'templates', 'models', 'packaged-models']
    project_metadata_filepath = 'project-data.yaml'

    now = datetime.datetime.utcnow()
    datestamp = now.strftime(MSMSeeder.core.datestamp_format_string)

    argparser = argparse.ArgumentParser(description='Initialize MSMSeeder project within the current directory. Creates necessary subdirectories and a project metadata .yaml file.')
    argparser.parse_args()

    # =========
    # Create necessary project directories
    # =========

    os.chdir(project_toplevel_dir)

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

    project_metadata = MSMSeeder.core.ProjectMetadata()
    init_metadata = {'init' : {'datestamp' : datestamp}}
    project_metadata.add_metadata(init_metadata)
    project_metadata.write(ofilepath=project_metadata_filepath)

    print 'Done.'

if __name__ == '__main__':
    init_msmseeder_project('.')

