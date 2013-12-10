MSMSeeder
=========

Generation of diverse structural ensembles from homology data, for the initialization of molecular dynamics simulation trajectories and construction of Markov state models.

Authors
-------

* John D. Chodera | john.chodera@choderalab.org
* Daniel L. Parton | danny.parton@choderalab.org
* Patrick B. Grinaway | patrick.grinaway@choderalab.org

Manifest
--------

targets/ - target sequences

templates/ - template sequences and structures

models/ - comparative models for initialization of simulations

simulations/ - simulation data

pylib/ - useful Python libraries

spark-distributed/ - Library with map/reduce-able functions

Dependencies
------------

* OpenMM - https://simtk.org/home/openmm
* Modeller - http://salilab.org/modeller/
* mpi4py - http://mpi4py.scipy.org/
* MDAnalysis - https://code.google.com/p/mdanalysis/
* Clustal Omega - http://www.clustal.org/omega/
* PyMOL (optional, for model alignment/visualization) - http://www.pymol.org/
* Stride (optional, for model alignment/visualization) - http://webclu.bio.wzw.tum.de/stride/ (also included as part of the VMD package - http://www.ks.uiuc.edu/Research/vmd/)
* Many other Python packages. Recommended aproach is to install either Enthought Canopy (https://www.enthought.com/products/canopy/) or Continuum Anaconda (https://store.continuum.io/)
* Apache Spark - http://spark.incubator.apache.org/
