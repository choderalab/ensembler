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
  ensembler refine_explicit
  ensembler package_models

The optional ``ensembler validate`` subcommand uses the
`MolProbity <http://molprobity.biochem.duke.edu/>`_ command-line tools to
conduct model quality validation based on criteria such as Ramachandran angles,
backbone distortion, and atom clashes.

The ``ensembler quickmodel`` subcommand allows the entire modeling pipeline to
be run in one go for a single target and a small number of templates. Note that
this command will not work with MPI.

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

Custom settings
===============

Many aspects of the behavior of Ensember can be specified by using the Python
API instead of the main command-line interface.  For API documentation, see the
source code, or view the docstrings in iPython.

Custom settings via the manual_overrides.yaml file
--------------------------------------------------

Some options can instead be specified via the ``manual_overrides.yaml`` file,
which is created when initializing a new Ensembler project. The file contains
an example configuration, with each line commented out. The user can thus
uncomment the relevant lines and edit as necessary.

::

  target-selection:
      domain-spans:
        ABL1_HUMAN_D0: 242-513
  template-selection:
      min-domain-len: 0
      max-domain-len: 350
      domain-spans:
          ABL1_HUMAN_D0: 242-513
      skip-pdbs:
          - 4CYJ
          - 4P41
          - 4P2W
          - 4QTD
          - 4Q2A
          - 4CTB
          - 4QOX
  refinement:
      ph: 8.0
      custom_residue_variants:
          DDR1_HUMAN_D0_PROTONATED:
              # keyed by 0-based residue index
              35: ASH

The above configuration makes the following specifications (in order of appearance):

 - Specifies a custom residue span for the target ``ABL1_HUMAN_D0``. This is useful in cases where a different domain span is desired from that annotated in UniProt.
 - Specifies minimum and maximum domain lengths for templates. Any domain with more than 350 residues would be excluded. The same custom residue span used for target domains is also specified for the template domains.
 - Certain PDB files can be skipped if they cause problems.
 - A custom pH level (default: 7.0) is set, which determines how protonation states are assigned by OpenMM prior to molecular dynamics refinement.
 - Custom residue variants are specified. This can be used to set specific protonation states, rather than rely purely on a defined pH level. These specified protonation states would override those determined by pH. The naming of residue variants (e.g. ``ASH``) follows the OpenMM conventions.

Additional Tools
================

Ensembler includes a ``tools`` submodule, which allows the user to conduct
various useful tasks which are not considered core pipeline functions. The
use-cases for many of these tools are quite specific, so they may not be
applicable to every project, and should be used with caution.

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

Generating unrefined model structures
-------------------------------------

In some cases it may be useful to analyze model structures which have not
undergone refinement, but which have topologies equivalent to the final refined
models. These structures are not saved by the main pipeline functions by
default, but can be regenerated using
``ensembler.tools.mktraj.MkTrajImplicitStart``. This code simply loads each
model structure with ``openmm``, adds hydrogens, and writes the resultant
structure as a pdb file (``implicit-start.pdb.gz``). It also combines the
structures into a trajectory (``traj-implicit-start.xtc``). This function is
accessed via the Python API as follows:

::

  from ensembler.tools.mktraj import MkTrajImplicitStart
  MkTrajImplicitStart(targetid='EGFR_HUMAN_D0')

