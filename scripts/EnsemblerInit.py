#!/usr/bin/env python
#
# Initialize project directory by creating necessary subdirectories and a project metadata file.
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 11 Mar 2014

import argparse
import ensembler.initproject


def main():
    # ========
    # Parse command-line arguments
    # ========

    argparser = argparse.ArgumentParser(
        description='Initialize Ensembler project by creating necessary subdirectories and a project metadata .yaml file.'
    )
    argparser.add_argument(
        '--project_dir', type=str, default='.',
        help='(Default: ".") Optionally provide a directory path in which to initialize the project.'
    )
    args = argparser.parse_args()

    ensembler.initproject.initproject(args.project_dir)


if __name__ == '__main__':
    main()