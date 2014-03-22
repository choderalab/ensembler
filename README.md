MSMSeeder
=========

Generation of diverse protein structural ensembles, for the initialization of molecular dynamics simulations and subsequent construction of Markov state models.

Authors
-------

* Daniel L. Parton | daniel.parton@choderalab.org
* John D. Chodera | john.chodera@choderalab.org
* Patrick B. Grinaway | patrick.grinaway@choderalab.org

Manifest
--------

scripts/ - Python wrapper scripts which accept parameters from the command-line or project metadata file

MSMSeeder/ - main code

tests/ - to test everything is working correctly, run nosetests from this top-level directory

Dependencies
------------

* OpenMM - https://simtk.org/home/openmm
* Modeller - http://salilab.org/modeller/
* mpi4py - http://mpi4py.scipy.org/
* mdtraj - http://mdtraj.org/
* PyMOL (optional, for model alignment/visualization) - http://www.pymol.org/
* Various other Python packages commonly used in scientific computing. Recommended aproach is to install either Enthought Canopy (https://www.enthought.com/products/canopy/) or Continuum Anaconda (https://store.continuum.io/)

