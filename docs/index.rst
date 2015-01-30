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

Ensembler maps a set of target sequences onto a set of template structures, and shepherds the resulting models through a series of successive refinement and filtering stages. This:

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

1. Retrieval of protein target sequences and template structures, e.g. from `UniProt <http://www.uniprot.org/>`_ and the `PDB <http://www.pdb.org>`_
2. (*optional*) Reconstruction of missing template loops
3. Model generation - each target sequence is mapped onto each available template structure
4. Culling of models based on close structural similarity
5. Refinement with implicit solvent molecular dynamics simulation
6. Solvation of models with explicit water
7. (*optional*) Refinement with explicit solvent molecular dynamics simulation
8. (*optional*) Packaging of models, ready for transfer or set-up on production simulation platforms such as ``Folding@Home``


License
-------
Ensembler is licensed under the GNU General Public License (GPL) v2.0.

.. raw:: html

   <div style="display:none">

.. toctree::
   :maxdepth: 2

   installation
   examples
   cli_usage

.. raw:: html

   </div>

