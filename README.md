Ensembler
=========

[![Binstar Badge](https://binstar.org/omnia/ensembler/badges/version.svg)](https://binstar.org/omnia/ensembler)
[![Documentation Status](https://readthedocs.org/projects/ensembler/badge/?version=latest)](http://ensembler.readthedocs.org/en/latest/)
[![Build Status](https://travis-ci.org/choderalab/ensembler.svg)](https://travis-ci.org/choderalab/ensembler)

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
2. Build models by mapping each target sequence onto every available template structure, using [Modeller](http://salilab.org/modeller/).
3. Filter out non-unique models (based on a RMSD cutoff).
4. Refine models with implicit solvent molecular dynamics simulation.
5. Refine models with explicit solvent molecular dynamics simulation.
6. (_optional_) Package and/or compress the final models, ready for transfer or for set-up on other platforms such as [Folding@Home](http://folding.stanford.edu/).

Installation
------------

Using conda (installs all dependencies except Modeller and Rosetta):

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
* mdtraj - http://mdtraj.org/
* MSMBuilder - http://msmbuilder.org/
* PDBFixer - https://github.com/pandegroup/pdbfixer
* BioPython
* NumPy
* lxml
* PyYAML
* docopt
* mock
* Optional:
  * Rosetta (optional, for template loop reconstruction) - https://www.rosettacommons.org/software
  * MPI4Py (allows many Ensembler functions to be run in parallel using MPI)
  * Pandas (required for certain analysis functions)
  * subprocess32 (if using Python 2)
  * PyMOL (optional, for model alignment/visualization) - http://www.pymol.org/

Recommended approach is to install using conda (https://store.continuum.io/cshop/anaconda/). This will install all dependencies except for Modeller, which must be installed separately by the user.
