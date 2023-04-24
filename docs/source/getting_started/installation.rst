============
Installation
============

Clone the repository and change directory:

.. code-block:: console

  $ git clone https://github.com/pnnl/isicle.git
  $ cd isicle/

Use `conda <https://www.anaconda.com/download/>`_ (or, more efficiently, `mamba <https://mamba.readthedocs.io/en/latest/>`_) to create a virtual environment with required dependencies.
Only 64-bit Linux and Mac (excluding `mobcal-shm <https://github.com/pnnl/mobcal-shm>`_) platforms are supported due to dependency requirements.

On Mac:

.. code-block:: console
  
  $ conda env create -f envs/osx.yml
  or
  $ mamba env create -f envs/osx.yml

On Linux:

.. code-block:: console
  
  $ conda env create -f envs/linux.yml
  or
  $ mamba env create -f envs/linux.yml

Activate the virtual environment:

.. code-block:: console
  
  $ conda activate deimos

Install ISiCLE using `pip <https://pypi.org/project/pip/>`_:

.. code-block:: console
  
  $ pip install -e .
