Download an executable for Linux or Mac via Conda
-------------------------------------------------

Pre-compiled LAMMPS binaries are available for macOS and Linux via the
`Conda <conda_>`_ package management system.

First, one must set up the Conda package manager on your system.  Follow the
instructions to install `Miniconda <mini_conda_install_>`_, then create a conda
environment (named `my-lammps-env` or whatever you prefer) for your LAMMPS
install:

.. code-block:: bash

   conda config --add channels conda-forge
   conda create -n my-lammps-env

Then, you can install LAMMPS on your system with the following command:

.. code-block:: bash

   conda activate my-lammps-env
   conda install lammps

The LAMMPS binary is built with the :ref:`KIM package <kim>`, which
results in Conda also installing the `kim-api` binaries when LAMMPS is
installed.  In order to use potentials from `openkim.org <openkim_>`_,
you can install the `openkim-models` package

.. code-block:: bash

   conda install openkim-models

If you have problems with the installation, you can post issues to `this
link <conda_forge_lammps_>`_.  Thanks to Jan Janssen
(Max-Planck-Institut fuer Eisenforschung) for setting up the Conda
capability.

.. _conda_forge_lammps: https://github.com/conda-forge/lammps-feedstock/issues
.. _openkim: https://openkim.org
.. _conda: https://docs.conda.io/en/latest/index.html
.. _mini_conda_install: https://docs.conda.io/en/latest/miniconda.html

.. note::

   If you have questions about these pre-compiled LAMMPS executables,
   you need to contact the people preparing those packages.  The LAMMPS
   developers have no control over their choices of how they configure
   and build their packages and when they update them.
