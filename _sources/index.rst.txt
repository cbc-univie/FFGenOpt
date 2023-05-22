.. FFGenOpt documentation master file, created by
   sphinx-quickstart on Tue May  9 09:23:55 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FFGenOpt's documentation!
====================================
**FFGenOpt** is a genetic algorithm for automatic MD forcefield parameter
refinement implemented in python. It targets quantum mechanical normal modes
and frequencies and tunes the input classical molecular dynamics force field
so that the normal modes in the MM and QM set match each other.

To get a quick rundown on the refinement procedure, you can take a look at
the Tutorial, explaining the process for generating parameters for
trifluoroacetate.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   documentation_pages/qm
   documentation_pages/forcefield
   documentation_pages/ffgenopt



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
