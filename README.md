MSMSeeder
=========

Software pipeline used to generate diverse protein structural ensembles, for
the purpose of seeding multiple parallel molecular dynamics simulations, and
subsequent construction of Markov state models. This branch contains a version of MSMSeeder which
runs via the distributed computing framework Spark (https://spark.apache.org/)

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

MSMSeeder/ - main code; can be used as a standalone Python library package

tests/ - to test whether the code is working correctly, run nosetests from this top-level directory

Installation
------------

    git clone https://github.com/choderalab/msmseeder.git
    cd msmseeder
    python setup.py install

Dependencies
------------

* OpenMM - https://simtk.org/home/openmm
* Modeller - http://salilab.org/modeller/
* mpi4py - http://mpi4py.scipy.org/
* mdtraj - http://mdtraj.org/
* Spark  - https://spark.apache.org/
* PyMOL (optional, for model alignment/visualization) - http://www.pymol.org/
* Various other Python packages commonly used in scientific computing. Recommended aproach is to install either Continuum Anaconda (https://store.continuum.io/) or Enthought Canopy (https://www.enthought.com/products/canopy/)


