.. _cli_usage:

**********************
Command-Line Interface
**********************

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

.. TODO note that you can print help with -h flag


init
====

::

  $ ensembler init

This sets up an Ensembler project in the current working directory. It creates
a number of directories and a metadata file (meta0.yaml).

gather_templates
================

::

  $ ensembler gather_targets

.. TODO ideally would generate this from docstrings
.. For now, could just give a brief overview of each command, and suggest people use the "-h" flag for further details
