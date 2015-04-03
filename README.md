Ensembler
=========

[![Binstar Badge](https://binstar.org/omnia/ensembler/badges/version.svg)](https://binstar.org/omnia/ensembler)
[![Documentation Status](https://readthedocs.org/projects/ensembler/badge/?version=latest)](http://ensembler.readthedocs.org/en/latest/)

Software pipeline for automating omics-scale protein modeling and simulation setup.

Online documentation
--------------------
Go to the [official online documentation](http://ensembler.readthedocs.org/).

Authors
-------

* Daniel L. Parton | daniel.parton@choderalab.org
* John D. Chodera | john.chodera@choderalab.org
* Patrick B. Grinaway | patrick.grinaway@choderalab.org

Overview of pipeline
--------------------

1. Retrieve protein target sequences and template structures.
2. Build models by mapping each target sequence onto every available template structure, using Modeller (http://salilab.org/modeller/).
3. Filter out non-unique models (based on a RMSD cutoff).
4. Refine models with implicit solvent molecular dynamics simulation.
5. Refine models with explicit solvent molecular dynamics simulation.
6. (_optional_) Package and/or compress the final models, ready for transfer or for set-up on other platforms such as [Folding@Home](http://folding.stanford.edu/).

Manifest
--------

scripts/ - Python wrapper scripts which accept parameters from the command-line or project metadata file

ensembler/ - main code; can be used as a standalone Python library package

tests/ - to test whether the code is working correctly, run nosetests from this top-level directory

Installation
------------

Using conda:

    conda config --add channels http://conda.binstar.org/omnia
    conda install ensembler

From source:

    git clone https://github.com/choderalab/ensembler.git
    cd ensembler
    python setup.py install

Dependencies
------------

* OpenMM - https://simtk.org/home/openmm
* Modeller - http://salilab.org/modeller/
* mpi4py - http://mpi4py.scipy.org/
* mdtraj - http://mdtraj.org/
* MSMBuilder - http://msmbuilder.org/
* PyMOL (optional, for model alignment/visualization) - http://www.pymol.org/
* Other Python packages: BioPython, mpi4py, numpy, lxml, PyYAML

Recommended approach is to install using conda (https://store.continuum.io/cshop/anaconda/). This will install all dependencies except for Modeller, which must be installed separately by the user.

Basic Usage
-----------

The package can be used either by running the command-line scripts listed
below, or by writing scripts to interact with the Ensembler Python library.
Each script can be run from the command-line with a '-h' flag, which will print
information on their usage. They are intended to be run in the following order:

1. EnsemblerInit
2. EnsemblerGatherTargets
3. EnsemblerGatherTemplates
5. EnsemblerRefineImplicitMD
6. EnsemblerSolvate
7. EnsemblerRefineExplicitMD
8. EnsemblerPackageModels

The scripts have been parallelized with MPI where useful. Furthermore, the MD
refinement steps use the OpenMM simulation package, which has particularly high
performance on GPUs.

Some scripts will write certain parameters to a project metadata file (named
"project-data.yaml" and stored in the project top-level directory). The purpose
of this is partly to keep a record of the input parameters used, but if wanted,
this file can be used to provide the same input parameters the next time the
same script is run. Note that input parameters read from command-line arguments
take precedence over those read from the project metadata file. A template
project metadata file is provided in the source code.

Sometimes it is necessary to manually override certain automatic aspects of the
target and template retrieval scripts. This can be achieved using a file named
'manual-specifications.yaml', to be placed in the project top-level directory
(a template is provided in the source code). Currently, this allows one to
manually specify individual domain spans for targets and templates retrieved
from UniProt, as well as minimum and maximum acceptable domain lengths.

Example commands
----------------

    EnsemblerInit
    EnsemblerGatherTargets --gather_from TargetExplorerDB --db_path ~/kinomeDB/database/database.xml
    EnsemblerGatherTemplates --gather_from UniProt --uniprot_query 'domain:"Protein kinase" AND reviewed:yes' --uniprot_domain_regex '^Protein kinase(?!; truncated)(?!; inactive)' --structure_paths ~/pdbfiles ~/siftsfiles
    EnsemblerBuildModels --targets ABL1_HUMAN_D0
    EnsemblerRefineImplicitMD --targets ABL1_HUMAN_D0
    EnsemblerSolvate --targets ABL1_HUMAN_D0
    EnsemblerRefineExplicitMD --targets ABL1_HUMAN_D0
    EnsemblerPackageModels --targets ABL1_HUMAN_D0 --package_for FAH

Additional notes on scripts
----------------

(run each script with '-h' flag for general info)

### EnsemblerGatherTargets and EnsemblerGatherTemplates

These scripts are used to provide target data (ids and sequences) and
template data (ids, sequences and structures). Two different
methods are provided for each script, specified with a command-line flag:

* _--gather\_from TargetExplorerDB_
    * Retrieves data from an existing
[TargetExplorer](https://github.com/choderalab/targetexplorerdb) database
* _--gather\_from UniProt_
    * retrieves data from the UniProt web server, as specified by two additional flags:
    * _--uniprot\_query_
        * Specifies the search string used to select UniProt entries
        * This uses the same syntax as the search function on the
UniProt site; not all types of syntax may be supported, but most basic searches
will work.
        * Note that this flag used on its own will return UniProt _entries_ (i.e. full-length proteins)
    * _--uniprot\_domain\_regex_
        * Pass a regular expression to this flag to subselect only certain target/template domains

It should also be noted that the output from these two scripts is simply a list
of target sequences (targets/targets.fa), a list of template sequences
(templates/templates.fa), and a set of template structures
(templates/structures/[id].pdb). The subsequent model-building and refinement
scripts can thus be run using any set of data specified in this way, allowing
the EnsemblerGatherTargets and EnsemblerGatherTemplates scripts to be bypassed or modified
if necessary. The main restrictions are that modified residues (e.g.
phosphotyrosines) cannot be included, and that ids and residues in the
templates.fa file must match with the template structure files.

