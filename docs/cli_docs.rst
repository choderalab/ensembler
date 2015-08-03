.. _cli_docs

************************************
Command-Line Interface Documentation
************************************

Overview
========

Ensembler can be used via the command-line tool ``ensembler`` or via the API.

For API documentation, see the source code.

The ``ensembler`` tool is operated via a number of subcommands, which should be executed successively ::

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

If desired, target-selection and template-selection can be set up manually,
rather than using the ``gather_targets`` and ``gather_templates`` subcommands.
Targets should be provided as a fasta-format file (``targets/targets.fa``)
containing target sequences and arbitrary identifiers.  Template sequences and
arbitrary identifiers should be provided in a fasta-format file
(``templates/templates-resolved-seq.fa``), and structures should be provided as
PDB-format coordinate files in the directory ``templates/structures-resolved``.
Each structure should be named XXX.pdb, where XXX matches the identifier in the
fasta file. The residues in the coordinate files should also match the
sequences in the fasta file.

Additional Tools
================

Ensembler includes a ``tools`` submodule, which allows the user to conduct
various useful tasks which are not considered core pipeline functions. The
use-cases for many of these tools are quite specific, so they may not be
applicable to every project, and should also be used with caution.

Residue renumbering according to UniProt sequence coordinates
-------------------------------------------------------------

::

  $ ensembler renumber_residues --target EGFR_HUMAN_D0

The given target ID must begin with a UniProt mnemonic, e.g. "EGFR_HUMAN".
This will output two files in the ``models/[target_id]`` directory:
``topol-renumbered-implicit.pdb`` and ``topol-renumbered-explicit.pdb``.
The coordinates are simply copied from the first example found for each of
``refined-implicit.pdb.gz`` and ``refined-explicit.pdb.gz``. The residue
numbers are renumbered according to the canonical isoform sequence coordinates
in the UniProt entry.
