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

For further details on their usage, see the main `command-line interface documentation <cli_docs.html>`_.

Example using the quickmodel function
=====================================

::

  $ ensembler quickmodel --target_uniprot_entry_name EGFR_HUMAN --uniprot_domain_regex '^Protein kinase' --template_pdbids 4KB8,4AF3 --no-loopmodel

Models human EGFR onto two templates selected via PDB IDs. The ``quickmodel`` function executes the entire modeling pipeline in one go, and is designed to work with only a few targets and templates. For generating larger numbers of models (such as entire protein families), the main pipeline functions should be used.

Example using the main pipeline functions
=========================================

::

  $ ensembler init

This sets up an Ensembler project in the current working directory. It creates
a number of directories and a metadata file (meta0.yaml).

::

  $ ensembler gather_targets --gather_from uniprot --query 'domain:"Protein kinase" AND taxonomy:9606 AND reviewed:yes' --uniprot_domain_regex '^Protein kinase(?!; truncated)(?!; inactive)'

Queries UniProt for all human protein kinases, and selects the domains of interest, as specified by the `regular expression <https://docs.python.org/2/library/re.html#regular-expression-syntax>`_ ("regex") passed to the final flag. At the time this documentation was written, five types of protein kinase domain were returned by the UniProt search, annotated as "Protein kinase", "Protein kinase; 1", "Protein kinase; 2", "Protein kinase; truncated", and "Protein kinase; inactive". The above regex selects the first three types of domain, and excludes the latter two. Sequences are written to a fasta file: ``targets/targets.fa``.

Targets are given IDs of the form ``[UniProt mnemonic]_D[domain id]``, which consists of the UniProt name for the target and an identifier for the domain (since a single target protein may contain multiple domains of interest). Example: ``EGFR_HUMAN_D0``.

::

  $ ensembler gather_templates --gather_from uniprot --query 'domain:"Protein kinase" AND reviewed:yes' --uniprot_domain_regex '^Protein kinase(?!; truncated)(?!; inactive)'

Queries UniProt for all protein kinases (of any species), selects the relevant domains, and retrieves sequence data and a list of associated PDB structures (X-ray and NMR only), which are then downloaded from the PDB. Template sequences are written in two forms - the first contains only residues resolved in the structure (``templates/templates-resolved-seq.fa``); the second contains the complete UniProt sequence containined within the span of the structure, including unresolved residues (``templates/templates-full-seq.fa``). Template structures (containing only resolved residues) are extracted and written to the directory ``templates/structures-resolved``. Templates containing the full sequences can optionally be generated with a subsequent step - the ``loopmodel`` function.

Templates are given IDs of the form ``[UniProt mnemonic]_D[domain id]_[PDB id]_[chain id]``, where the final two elements represent the PBD ID and chain identifier. Example: ``EGFR_HUMAN_D0_2GS7_B``.

::

  $ ensembler loopmodel

(*Optional*)
Reconstruct template loops which were not resolved in the original PDB structure, using ``Rosetta loopmodel``. This tends to result in higher quality models. The reconstructed template structures are written to the directory ``templates/structures-modeled-loops``.

::

  $ ensembler align

Conducts pairwise alignments of target sequences against template sequences. These alignments are used to guide the subsequent modeling step, and are stored in directories of the form ``models/[target ID]/[template ID]/alignment.pir``. The ``.pir`` alignment format is an ascii-based format required by ``Modeller``.

If the ``loopmodel`` function was used previously, then templates which have been successfully remodeled will be selected for this alignment and the subsequent modeling steps. Otherwise, Ensembler defaults to using the template structures which contain only resolved residues.

::

  $ ensembler build_models

Creates models by mapping each target sequence onto each template structure, using the ``Modeller`` `automodel <https://salilab.org/modeller/manual/node15.html>`_ function.

::

  $ ensembler cluster
