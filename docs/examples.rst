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

Retrieves sequence data for all human protein kinases listed in UniProt, and selects only domains labeled "Protein kinase", "Protein kinase; 1" or "Protein kinase; 2" (not "Protein kinase; truncated" or "Protein kinase; inactive"). Sequences are written to a fasta file: `targets/targets.fa`. Targets are given IDs of the form `[UniProt mnemonic]_D[domain id]`, which consists of the UniProt name for the target and an identifier for the domain (since a single target protein may contain multiple domains of interest). Example: `EGFR_HUMAN_D0`.

::

  $ ensembler gather_templates --gather_from uniprot --query 'domain:"Protein kinase" AND reviewed:yes' --uniprot_domain_regex '^Protein kinase(?!; truncated)(?!; inactive)'

Queries UniProt for all protein kinases (of any species), selects the relevant domains, and retrieves sequence data and a list of associated PDB structures (X-ray and NMR only), which are then downloaded from the PDB. Finally, template structures are extracted and written to the directory `templates/structures-resolved`.

::

  $ ensembler loopmodel

(*Optional*)
Uses `Rosetta loopmodel` to reconstruct template loops which were not resolved in the original PDB structure. This tends to result in higher-quality models, following the subsequent build_models step. The reconstructed template structures are written to the directory `templates/structures-modeled-loops`.

::

  $ ensembler build_models

Uses `Modeller` to map each target sequence onto each template structure.
