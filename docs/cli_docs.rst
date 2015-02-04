.. _cli_docs

************************************
Command-Line Interface Documentation
************************************

Overview
========

Ensembler can be used via the command-line tool ``ensembler`` or via the API.

For API documentation, see the source code.

The ``ensembler`` tool is operated via a number of subcommands: ::

  ensembler init
  ensembler gather_targets
  ensembler gather_templates
  ensembler loopmodel
  ensembler align
  ensembler build_models
  ensembler cluster
  ensembler refine_implicit
  ensembler solvate
  ensembler refine_implicit
  ensembler package_models

Furthermore, the ``ensembler quickmodel`` subcommand allows the entire modeling
pipeline to be run in one go for a single target and a small number of
templates. Note that this command will not work with MPI.

To print helpstrings for each subcommand, pass the ``-h`` flag.
