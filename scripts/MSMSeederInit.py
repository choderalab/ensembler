#!/usr/bin/env python
#
# Initialize project directory by creating necessary subdirectories and a project metadata file.
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 11 Mar 2014

import os
import argparse
import msmseeder.initproject

# ========
# Parse command-line arguments
# ========

def parse_args():
    argparser = argparse.ArgumentParser(description='Initialize MSMSeeder project by creating necessary subdirectories and a project metadata .yaml file.')
    argparser.add_argument('--project_dir', type=str, default='.', help='(Default: ".") Optionally provide a directory path in which to initialize the project.')
    args = argparser.parse_args()
    return args

def main(args):
    project_dir = os.path.abspath(args.project_dir)

    if not os.path.exists(project_dir):
        os.makedirs(project_dir)

    msmseeder.initproject.init(project_dir)

if __name__ == '__main__':
    args = parse_args()
    main(args)
