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

`conda <http://www.continuum.io/blog/conda>`_ is a package manager built for scientific Python. Unlike ``easy_install`` or ``pip``, it handles binaries and binary dependencies, which are critical for most scientific workflows.

``conda`` can be obtained by installing `Continuum Anaconda <https://store.continuum.io/>`_ - an awesome free Python distribution for scientific computing. The standard installation contains many of Ensembler's dependencies. Alternatively, for a more minimal Python set-up you can install `miniconda <http://conda.pydata.org/miniconda.html>`_, which contains only ``conda`` and Python.

To install Ensembler with ``conda``, use the following commands ::

  $ conda config --add channels http://conda.binstar.org/dannyparton
  $ conda install ensembler

``conda`` will automatically install all dependencies except for `Modeller <https://salilab.org/modeller/>`_ and `Rosetta <https://www.rosettacommons.org/software>`_ (optional). These require licenses (free for academic non-profit use), and will have to be installed according to the instructions for those packages.

.. warning:: If you are running Ensembler on OS X and see errors relating to OpenMPI, e.g. ``Sorry!  You were supposed to get help about: opal_init:startup:internal-failure``), this is a `known issue with Anaconda <https://github.com/ContinuumIO/anaconda-issues/issues/96>`_. A simple workaround is to create a symlink from ``/opt/anaconda1anaconda2anaconda3`` to your Anaconda installation, e.g. ``sudo ln -s ~/anaconda /opt/anaconda1anaconda2anaconda3``.


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

    `MSMBuilder <http://msmbuilder.org/>`_
        Statistical models for biomolecular dynamics.

    `PDBFixer <https://github.com/pandegroup/pdbfixer>`_
        PDB structure modeling

    `BioPython <http://biopython.org/wiki/Main_Page>`_
        Collection of Python tools for computational biology and
        bioinformatics.

    `NumPy <http://numpy.scipy.org/>`_
        Numpy is the base package for numerical computing in Python.

    `lxml <http://lxml.de/>`_
        For working with XML files.

    `PyYAML <http://pyyaml.org/>`_
        For working with YAML files.

    `docopt <http://docopt.org/>`_
        For building command-line interfaces.

    `mock <http://www.voidspace.org.uk/python/mock/>`_
        For testing in Python

Optional packages:

    `MPI4Py <http://mpi4py.scipy.org/>`_
        Allows many Ensembler functions to be run in parallel using MPI.

    `Rosetta <https://www.rosettacommons.org/software>`_
        Protein modeling suite. The ``loopmodel`` function is optionally used
        by Ensembler to reconstruct missing loops in template structures.

    `subprocess32 <https://pypi.python.org/pypi/subprocess32/>`_
        (If running Python 2.)
        Backport of the Python 3 subprocess module for Python 2. Used to run
        command-line programs such as Rosetta loopmodel. Includes timeout
        functionality which is particularly useful for the template loop
        reconstruction routine.

    `Pandas <http://pandas.pydata.org>`_
        Some functionality, including the ``quickmodel`` and ``inspect``
        functions, requires pandas.

Manually Installing the Dependencies
------------------------------------

Linux
++++++
If you're on ubuntu and have root, you can install most dependencies through your package manager (``apt-get``). ::

  $ sudo apt-get install python-dev

Mac
+++
If you're on mac and want a package manager, you should be using `homebrew <http://mxcl.github.io/homebrew/>`_ and ``brews``'s Python (see `this page <https://github.com/mxcl/homebrew/wiki/Homebrew-and-Python>`_ for details). For example, numpy can be installed with ``brew`` as follows: ::

  $ brew tap Homebrew/python
  $ brew install python
  $ brew install numpy

Then, you can install many of the remaining packages with ``pip``. ::

  $ pip install lxml

Windows
+++++++
Chris Gohlke maintains windows binary distributions for an ever-growing
set of Python extensions on `his website <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_.
Download and install the the installers for setuptools, nose, numpy, scipy, numexpr, pandas and tables.

Testing Your Installation
=========================
Running the tests is a great way to verify that everything is working. The test
suite uses `nose <https://nose.readthedocs.org/en/latest/>`_ and `mock
<http://www.voidspace.org.uk/python/mock/>`_, which you can pick up via
``conda`` or ``pip`` if you don't already have them. ::

  $ conda install nose mock

To run the unit tests: ::

  $ nosetests ensembler -a unit

Further tests are available which check interoperation of Ensembler with
software dependencies such as Modeller and Rosetta loopmodel, or with external public
databases such as UniProt, or are excluded from the unit tests due to being
slow. To run them: ::

  $ nosetests ensembler -a non_conda_dependencies -a network -a slow
