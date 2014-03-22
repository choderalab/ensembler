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

MSMSeeder/ - main code; can be used as a standalone library package

tests/ - to test everything is working correctly, run nosetests from this top-level directory

Basic Usage
-----------

Make sure this directory is visible to your Python interpreter (e.g. by adding
it to your $PYTHONPATH environment variable).

You may also find it useful to add the scripts/ directory to your $PATH
environment variable, so you can execute the contained scripts directly.

The recommended approach for familiarizing yourself with this package is to
first try out the scripts, but the library functions can also be used directly
if wanted. The scripts can each be run from the command-line with a '-h' flag,
which will print information on their usage. They are intended to be run in the
following order:

1. InitMSMSeederProject.py
2. GatherTargets.py
3. GatherTemplates.py
4. BuildModels.py
5. RefineImplicitMD.py
6. Solvate.py
7. RefineExplicitMD.py

Further scripts can be used to package and compress the results, ready for
transfer or for set-up on the Folding@Home distributed computing platform.

Dependencies
------------

* OpenMM - https://simtk.org/home/openmm
* Modeller - http://salilab.org/modeller/
* mpi4py - http://mpi4py.scipy.org/
* mdtraj - http://mdtraj.org/
* PyMOL (optional, for model alignment/visualization) - http://www.pymol.org/
* Various other Python packages commonly used in scientific computing. Recommended aproach is to install either Enthought Canopy (https://www.enthought.com/products/canopy/) or Continuum Anaconda (https://store.continuum.io/)

