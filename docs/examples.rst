.. _examples:

**************
Usage Examples
**************

There are two main ways to use the Ensembler command-line interface. The
``quickmodel`` function performs the entire modeling pipeline in one go, and is
designed to work with only a few targets and templates. For generating larger
numbers of models (such as entire protein families), the main pipeline
functions should be used. These perform each stage of the modeling process
individually, and the most computationally intensive stages can be run in
parallel to increase performance.

For further details on their usage, see the main :ref:`command-line interface documentation <cli_docs>`.

Example using the quickmodel function
=====================================

::

  $ ensembler quickmodel --target_uniprot_entry_name EGFR_HUMAN --uniprot_domain_regex '^Protein kinase' --template_pdbids 4KB8,4AF3 --no-loopmodel

Models human EGFR onto two templates selected via PDB IDs.


Example using the main pipeline functions
=========================================

::

  $ ensembler init

This sets up an Ensembler project in the current working directory. It creates
a number of directories and a metadata file (meta0.yaml).

::

  $ ensembler gather_targets --gather_from uniprot --query 'domain:"Protein kinase" AND taxonomy:9606 AND reviewed:yes' --uniprot_domain_regex '^Protein kinase(?!; truncated)(?!; inactive)'

Gathers sequence data from UniProt for protein kinases, and selects only domains labeled "Protein kinase", "Protein kinase; 1" or "Protein kinase; 2" (not "Protein kinase; truncated" or "Protein kinase; inactive").
