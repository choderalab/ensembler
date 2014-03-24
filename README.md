MSMSeeder
=========

Software pipeline used to generate diverse protein structural ensembles, for
the purpose of seeding multiple parallel molecular dynamics simulations, and
subsequent construction of Markov state models.

Authors
-------

* Daniel L. Parton | daniel.parton@choderalab.org
* John D. Chodera | john.chodera@choderalab.org
* Patrick B. Grinaway | patrick.grinaway@choderalab.org

Overview
--------

The first stage in the software pipeline is the retrieval of protein target
sequences and protein template structures. Models are subsequently built by
mapping each target sequence onto every available template structure using
Modeller (http://salilab.org/modeller/). After filtering out non-unique models
(based on an RMSD cutoff), the models are subjected to a two-stage molecular
dynamics simulation refinement; the first with implicit solvation, the second
with explicit solvation. For the explicit solvation step, the models are
solvated such that each system contains the same number of waters. Finally,
scripts are provided to package and compress the resulting models, ready for
transfer or for set-up on other platforms such as Folding@Home.

Manifest
--------

scripts/ - Python wrapper scripts which accept parameters from the command-line or project metadata file

MSMSeeder/ - main code; can be used as a standalone Python library package

tests/ - to test whether the code is working correctly, run nosetests from this top-level directory

Installation
------------

Clone this repository, move into this directory, then run:

    python setup.py install

Dependencies
------------

* OpenMM - https://simtk.org/home/openmm
* Modeller - http://salilab.org/modeller/
* mpi4py - http://mpi4py.scipy.org/
* mdtraj - http://mdtraj.org/
* PyMOL (optional, for model alignment/visualization) - http://www.pymol.org/
* Various other Python packages commonly used in scientific computing. Recommended aproach is to install either Enthought Canopy (https://www.enthought.com/products/canopy/) or Continuum Anaconda (https://store.continuum.io/)

Basic Usage
-----------

The package can be used either by running the command-line scripts listed
below, or by writing scripts to interact with the MSMSeeder Python library. The
former approach is recommended for initially familiarizing yourself with the
software. Each script can be run from the command-line with a '-h' flag, which
will print information on their usage. They are intended to be run in the
following order:

1. InitMSMSeederProject.py
2. GatherTargets.py
3. GatherTemplates.py
4. BuildModels.py
5. RefineImplicitMD.py
6. Solvate.py
7. RefineExplicitMD.py
8. PackageModels.py

The scripts have been parallelized with MPI where useful. Furthermore, the MD
refinement steps use the OpenMM simulation package, which has particularly high
performance on GPUs.

Example commands
----------------

    InitMSMSeederProject.py
    GatherTargets.py --gather_from TargetExplorerDB --db_path ~/kinomeDB/database/database.xml
    GatherTemplates.py --gather_from UniProt --uniprot_query 'domain:"Protein kinase" AND reviewed:yes' --uniprot_domain_regex '^Protein kinase(?!; truncated)(?!; inactive)' --structure_paths ~/pdbfiles ~/siftsfiles
    BuildModels.py --targets ABL1_HUMAN_D0
    RefineImplicitMD.py --targets ABL1_HUMAN_D0
    Solvate.py --targets ABL1_HUMAN_D0
    RefineExplicitMD.py --targets ABL1_HUMAN_D0
    PackageModels.py --targets ABL1_HUMAN_D0 --package_for FAH

Additional notes on scripts
----------------

(run each script with '-h' flag for general info)

### GatherTargets.py and GatherTemplates.py

These two scripts are used to provide target data (ids and sequences) and
template data (ids, sequences and structures). In each case, two different
methods are provided, specified using the *--gather\_from* flag. The
'TargetExplorerDB' method allows one to retrieve data from an existing
TargetExplorer database (https://github.com/choderalab/targetexplorerdb), while
the 'UniProt' method retrieves data from the UniProt web server. In the latter
case, the *--uniprot\_query* flag specifies the search string used to select
UniProt entries. This uses the same syntax as the search function on the
UniProt site; not all types of syntax may be supported, but most basic searches
will work. While the initial data returned is in the form of UniProt *entries*
(i.e.  full-length proteins), the targets and templates selected by the scripts
must be protein domains. Without further specification, the scripts will select
all protein domains contained within the returned Uniprot entries. Often it
will be necessary to select only a subset of these protein domains - e.g. if
protein kinase domains are the desired target, then it will be necessary to
filter out other domains contained in the full-length kinase proteins - and
this can be achieved by supplying an appropriate regular expression along with
the *--uniprot\_domain\_regex* flag.

It should also be noted that the output from these two scripts is simply a list
of target sequences (targets/targets.fa), a list of template sequences
(templates/templates.fa), and a set of template structures
(templates/structures/[id].pdb). The subsequent model-building and refinement
scripts can thus be run using any set of data specified in this way, allowing
the GatherTargets.py and GatherTemplates.py scripts to be bypassed or modified
if necessary. The main restrictions are that modified residues (e.g.
phosphotyrosine) cannot be included, and that ids and residues in the
templates.fa file must match with the template structure files.

