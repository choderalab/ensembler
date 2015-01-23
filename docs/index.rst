.. Ensembler documentation master file, created by
   sphinx-quickstart on Fri Jan 16 17:34:52 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Ensembler
=====================================

*Software pipeline for generating diverse protein configurational ensembles
suitable for seeding highly parallel molecular simulations.*


Overview of Pipeline
--------------------

1. Gather protein target sequences and template structures.
2. (*optional*) Reconstruct missing template loops, using ``Rosetta loopmodel``.
3. Build models by mapping each target sequence onto each available template structure, using ``Modeller``.
4. Culling of models based on close structural similarity.
5. Refinement with implicit solvent molecular dynamics simulation, using ``OpenMM``.
6. Solvate models.
7. (*optional*) Refinement with explicit solvent molecular dynamics simulation, using ``OpenMM``.
8. (*optional*) Package models, ready for transfer or set-up on other platforms such as ``Folding@Home``.


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

