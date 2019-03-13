ISiCLE
======

Overview
--------
ISiCLE, or the _in silico_ chemical library engine, is a pipeline for high-accuracy chemical property calculation. ISiCLE takes an [InChI](https://en.wikipedia.org/wiki/International_Chemical_Identifier) or [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) string as input, generates an initial 3D conformation, and subsequently optimizes this initial structure through molecular dynamics simulations and quantum chemistry optimizations. Finally, ISiCLE simulates desired properties (e.g. collision cross section, NMR chemical shifts) for each conformer yielded during molecular dynamics simulations to produce a single value, Boltzmann-weighted by relative Gibb's free energy, giving emphasis to properties from highly probable conformations.

ISiCLE is implemented using the [Snakemake](https://snakemake.readthedocs.io) workflow management system, enabling scalability, portability, provenance, fault tolerance, and automatic job restarting. Snakemake provides a readable Python-based workflow definition language and execution environment that scales, without modification, from single-core workstations to compute clusters through as-available job queuing based on a task dependency graph.

<p align="center">
  <img align="center" src="resources/schematic.svg" width="40%" height="40%">
</p>

Installation
------------
Use [``conda``](https://www.anaconda.com/download/) to create a new virtual environment with required dependencies:
```bash
conda create -n isicle -c bioconda -c openbabel -c rdkit -c ambermd python=3.6.1 openbabel rdkit ambertools snakemake numpy pandas yaml statsmodels
```

Additionally, ensure the following third-party software is installed and added to your ``PATH``:
* [cxcalc](https://chemaxon.com/products/marvin/download) (license required)
* [NWChem](http://www.nwchem-sw.org/index.php/Download) (not required for CCS _Lite_)

Activate the virtual environment:
```
conda activate isicle
```

Install ISiCLE using [``pip``](https://pypi.org/project/pip/):
```bash
# clone/install
git clone https://github.com/pnnl/isicle.git
pip install isicle/

# direct
pip install git+https://github.com/pnnl/isicle
```

Getting Started
---------------
For usage overview, use ``isicle --help`` or ``-h``. Currently, available modules include ``prep`` for input preparation, ``ccs`` for collision cross section calculation, and ``shifts`` for NMR chemical shift calculation. For all modules, a Snakemake configuration file in [YAML](http://yaml.org/) format is required. ISiCLE will try to find ``config.yaml`` in the current directory, else a configuration file must be specified through the ``--config`` flag. Default [workflow](resources/example_config.yaml) and [cluster](resources/example_cluster.yaml) configurations are provided, but these are intended to be modified and supplied by the user to accomodate workflow-specific needs.

For the ``prep`` module, ISiCLE assumes the user starts with a text file with InChI or SMILES strings on each line. This ensures each input has a unique filname based on its InChI key identifier. We recommend using SMILES, as in some instances InChI processing can lead to unexpected results, though these occurences are rare. Detailed instructions can be accessed through the help flag (``isicle prep --help`` or ``-h``).
```bash
isicle prep input_list.txt
```

For the ``ccs`` module, the user must specify calculation mode (``lite`` or ``standard``), followed by any additional flags (see ``isicle ccs --help`` or ``-h``). Before beginning a simulation, we recommend use of the ``--dryrun`` flag to ensure the run is configured correctly. For desktop environments, we recommend using ``lite`` mode:
```bash
isicle ccs lite --cores 4 --dryrun
```

For ``slurm`` cluster environments, ``standard`` mode can be used:
```bash
isicle ccs standard --cluster cluster.yaml --jobs 999 --dryrun
```

The ``shifts`` module does not require selection of a calculation mode, but is otherwise configured the same way as the ``ccs`` module. See ``isicle shifts --help`` or ``-h`` for a full list of options. We recommend use of supercomputing resources for the ``shifts`` module:
```bash
isicle shifts --cluster cluster.yaml --jobs 999 --dryrun
```

Citing ISiCLE
-------------
If you would like to reference ISiCLE in an academic paper, we ask you include the following references:

* Colby, S.M., Thomas, D.G., Nunez, J.R., Baxter, D.J., Glaesemann, K.R., Brown, J.M., Pirrung, M.A., Govind, N., Teeguarden, J.G., Metz, T.O., and Renslow, R.S., 2018. _ISiCLE: A molecular collision cross section calculation pipeline for establishing large in silico reference libraries for compound identification_. arXiv preprint arXiv:1809.08378.
* Yesiltepe, Y., Nunez, J.R., Colby, S.M., Thomas, D.G., Borkum, M.I., Reardon, P.N., Washton, N.M., Metz, T.O., Teeguarden, J.T., Govind, N., and Renslow, R.S., 2018. _An automated framework for NMR chemical shift calculations of small organic molecules_. Journal of Cheminformatics. In press.
* ISiCLE, version 0.1.0 http://github.com/pnnl/isicle (accessed Oct 2018)

The first is a [preprint paper](https://arxiv.org/abs/1809.08378) describing ISiCLE for CCS, the second describes ISiCLE for NMR chemical shifts, and the third is to cite the software package (update version and access date appropriately).

Disclaimer
----------
This material was prepared as an account of work sponsored by an agency of the United States Government. Neither the United States Government nor the United States Department of Energy, nor Battelle, nor any of their employees, nor any jurisdiction or organization that has cooperated in the development of these materials, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness or any information, apparatus, product, software, or process disclosed, or represents that its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or any agency thereof, or Battelle Memorial Institute. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or any agency thereof.

PACIFIC NORTHWEST NATIONAL LABORATORY operated by BATTELLE for the UNITED STATES DEPARTMENT OF ENERGY under Contract DE-AC05-76RL01830
