.. _installation:

************
Installation
************

Platforms
=========

Ensembler is developed and tested on Linux and Mac. Windows is not well supported.

Supported Hardware
------------------
The simulation refinement stages of Ensembler use OpenMM, which performs best on CUDA-enabled GPUs. However, OpenMM also has an optimized CPU platform, so use of GPUs is optional.


Install with Conda
------------------
.. _install-with-conda:

`conda <http://www.continuum.io/blog/conda>`_ is a python package manager built for scientific python. Unlike ``easy_install`` or ``pip``, it handles binaries and binary dependencies, which are critical for most scientific workflows. If you're a ``conda`` user, you can install Ensembler by adding the relevant channel. If you're not a conda user, you should look into it.

To install Ensembler with conda, use the following commands ::

  $ conda config --add channels http://conda.binstar.org/dannyparton
  $ conda install ensembler

.. note:: ``conda`` will automatically install many of the tricky dependencies from binary packages automatically! The easiest way to get conda is with the `Anaconda python distribution <https://store.continuum.io/cshop/anaconda/>`_.


Install from Source
-------------------
Clone the source code repository from github ::

  $ git clone git://github.com/choderalab/ensembler.git

Then, in the directory containing the source code, you can install it with. ::

  $ python setup.py install

Dependencies
============

To use Ensembler, the following libraries and software will need to be installed.

    Linux, Mac OS X or Windows operating system
        We develop mainly on 64-bit Linux and Mac machines. Windows is not
        well supported.

    `Python <http://python.org>`_ >= 2.6
        The development package (``python-dev`` or ``python-devel``
        on most Linux distributions) is recommended.

    `Modeller <https://salilab.org/modeller/>`_
        Comparative modeling of protein structures.

    `OpenMM <https://simtk.org/home/openmm>`_
        Molecular simulation toolkit, with GPU-accelerated simulation platform.

    `MDTraj <http://mdtraj.org/>`_
        Simulation trajectory analysis library.

    `BioPython <http://biopython.org/wiki/Main_Page>`_
        Collection of Python tools for computational biology and
        bioinformatics.

    `NumPy <http://numpy.scipy.org/>`_
        Numpy is the base package for numerical computing in python.

    `lxml <http://lxml.de/>`_
        For working with XML files.

    `PyYAML <http://pyyaml.org/>`_
        For working with YAML files.

    `docopt <http://docopt.org/>`_
        For building command-line interfaces.

Optional packages:

    `MPI4Py <http://mpi4py.scipy.org/>`_
        Allows many Ensembler functions to be run in parallel using MPI.

    `Rosetta <https://www.rosettacommons.org/software>`_
        Protein modeling suite. The ``loopmodel`` function is optionally used
        by Ensembler to reconstruct missing loops in template structures.

    `subprocess32 <https://pypi.python.org/pypi/subprocess32/>`_
        Backport of the Python 3 subprocess module for Python 2. Used to run
        command-line programs such as Rosetta loopmodel. Includes timeout
        functionality which is particularly useful for the template loop
        reconstruction routine.

    `Pandas <http://pandas.pydata.org>`_
        Some functionality, including the ``quickmodel`` and ``inspect``
        functions, requires pandas.

Avoid Hassles with Anaconda or Canopy
-------------------------------------

An easy way to obtain many of the dependencies is to install one of the
pre-packaged scientific python distributes like `Enthought's Canopy
<https://www.enthought.com/products/canopy/>`_ or `Continuum's Anaconda
<https://store.continuum.io/>`_. These distributions already contain all of the
dependences, and are distributed via 1-click installers for Windows, Mac and
Linux.

.. note:: The developers personally recommend Continuum's Anaconda. It's free, includes the **AWESOME** conda package manager, and quite simple to use.

Modeller and Rosetta (optional) are the only dependencies not available via
these distributions, and will have to be installed according to the
instructions for those packages.

Manually Installing the Dependencies
------------------------------------

Linux
++++++
If you're on ubuntu and have root, you can install most dependencies through your package manager (``apt-get``). ::

  $ sudo apt-get install python-dev

Mac
+++
If you're on mac and want a package manager, you should be using `homebrew <http://mxcl.github.io/homebrew/>`_ and ``brews``'s python (see `this page <https://github.com/mxcl/homebrew/wiki/Homebrew-and-Python>`_ for details). For example, numpy can be installed with ``brew`` as follows: ::

  $ brew tap Homebrew/python
  $ brew install python
  $ brew install numpy

Then, you can install many of the remaining packages with ``pip``. ::

  $ pip install lxml

Windows
+++++++
Chris Gohlke maintains windows binary distributions for an ever-growing
set of python extensions on `his website <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_.
Download and install the the installers for setuptools, nose, numpy, scipy, numexpr, pandas and tables.

Testing Your Installation
=========================
Running the tests is a great way to verify that everything is working. The test
suite uses `nose <https://nose.readthedocs.org/en/latest/>`_, which you can
pick up via ``conda`` or ``pip`` if you don't already have it. ::

  $ conda install nose

Currently, the best way to run the tests is to go to the Ensembler installation
directory (e.g.
``~/anaconda/lib/python2.7/site-packages/ensembler-0.2-py2.7.egg/ensembler``) and
run the unit tests with: ::

  $ nosetests -a unit

There is also a suite of integration tests, which test interoperation of
Ensembler with software dependencies such as Modeller and Rosetta loopmodel, or
external databases such as UniProt. Note that many of these tests run much more slowly
than the unit tests. To run them: ::

  $ nosetests -a integration
