MSMSeeder
=========

Generation of diverse protein structural ensembles, for the initialization of molecular dynamics simulations and subsequent construction of Markov state models.

Authors
-------

* John D. Chodera | john.chodera@choderalab.org
* Daniel L. Parton | danny.parton@choderalab.org
* Patrick B. Grinaway | patrick.grinaway@choderalab.org

Manifest
--------

scripts/ - Python wrapper scripts which accept parameters from the command-line or project metadata file

MSMSeeder/ - main code

tests/ - to test everything is working correctly, run nosetests from this top-level directory

models/scripts/ - deprecated; old scripts in the process of being refactored into the new codebase

Dependencies
------------

* OpenMM - https://simtk.org/home/openmm
* Modeller - http://salilab.org/modeller/
* mpi4py - http://mpi4py.scipy.org/
* MDAnalysis - https://code.google.com/p/mdanalysis/
* PyMOL (optional, for model alignment/visualization) - http://www.pymol.org/
* Stride (optional, for model alignment/visualization) - http://webclu.bio.wzw.tum.de/stride/ (also included as part of the VMD package - http://www.ks.uiuc.edu/Research/vmd/)
* Many other Python packages. Recommended aproach is to install either Enthought Canopy (https://www.enthought.com/products/canopy/) or Continuum Anaconda (https://store.continuum.io/)

