.. Ensembler documentation master file, created by
   sphinx-quickstart on Fri Jan 16 17:34:52 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Ensembler
=====================================

.. *Ensembler generates diverse protein configurational ensembles suitable for seeding highly parallel molecular simulations.*


*Ensembler is a Python application and library which generates diverse protein
configurational ensembles suitable for seeding highly parallel molecular
simulations.*

Ensembler:
 - Automates the time-consuming process of setting up protein systems for molecular simulation
 - Exploits the entire variety of available genomic and structural data to provide diverse arrays of configurational models
 - Allows the study of proteins across entire (super)families

The generated models can be used as starting configurations for highly parallel
molecular simulations, which take advantage of the increasing availability of
parallel or distributed compute architectures. This approach is particularly
beneficial when used in conjunction with recent techniques for constructing
kinetic models, such as Markov state models, which can generate insight from
multiple independent simulation trajectories.

Ensembler can be run on a single computer or on a parallel compute cluster, and makes use of a number of external packages:
 - `Modeller <https://salilab.org/modeller>`_ for comparative modeling of target sequences onto template structures
 - `Rosetta <https://www.rosettacommons.org/software>`_ loopmodel for reconstruction of missing template loops
 - `OpenMM <https://simtk.org/home/openmm>`_ for model refinement with highly efficient, GPU-acclerated, molecular dynamics simulation
 - `MDTraj <http://mdtraj.org>`_ for ultra-fast RMSD calculation


Overview of the Ensembler Pipeline
----------------------------------

1. Gather protein target sequences and template structures
2. (*optional*) Reconstruct missing template loops, using ``Rosetta loopmodel``
3. Build models by mapping each target sequence onto each available template structure, using ``Modeller``
4. Culling of models based on close structural similarity
5. Refinement with implicit solvent molecular dynamics simulation, using ``OpenMM``
6. Solvate models
7. (*optional*) Refinement with explicit solvent molecular dynamics simulation, using ``OpenMM``
8. (*optional*) Package models, ready for transfer or set-up on other platforms such as ``Folding@Home``


License
-------
Ensembler is licensed under the GNU General Public License (GPL) v2.0.

.. toctree::
   :maxdepth: 2

   installation
   basic_example
   cli_usage

.. raw:: html

   <div style="display:none">

.. raw:: html

   </div>

