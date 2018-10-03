ISiCLE
======

Overview
--------
ISiCLE, or the _in silico_ chemical library engine, is a pipeline for high-accuracy chemical property calculation. ISiCLE takes an InChI (international chemical identifier) string as input, generates an initial 3D conformation, and subsequently optimizes this initial structure through molecular dynamics simulations and quantum chemistry optimizations. Finally, ISiCLE simulates desired properties (e.g. collision cross section, CCS) for each conformer yielded during molecular dynamics simulations to produce a single value, Boltzmann-weighted by Gibb's free energy, giving emphasis to properties from highly probable conformations.

ISiCLE is implemented using the snakemake workflow management system, enabling scalability, portability, provenance, fault tolerance, and automatic job restarting. Snakemake provides a readable Python-based workflow definition language and execution environment that scales, without modification, from single-core workstations to compute clusters through as-available job queuing based on a task dependency graph.

Installation
------------
Simply clone ISiCLE to your workstation or cluster, ensuring the following Python packages are installed:
* Snakemake
* OpenBabel, PyBel


Getting Started
---------------


Citing ISiCLE
-------------
If you would like to reference ISiCLE in an academic paper, we ask you use both of the following references::

* Colby, S.M., Thomas, D.G., Nunez, J.R., Baxter, D.J., Glaesemann, K.R., Brown, J.M., Pirrung, M.A., Govind, N., Teeguarden, J.G., Metz, T.O. and Renslow, R.S., 2018. ISiCLE: A molecular collision cross section calculation pipeline for establishing large in silico reference libraries for compound identification. arXiv preprint arXiv:1809.08378.
* ISiCLE, version 0.1.0 http://github.com/pnnl/isicle (accessed Oct 2018)

The first is a preprint paper describing ISiCLE, the second is to cite the software package (update version and access date appropriately).

Disclaimer
----------
This material was prepared as an account of work sponsored by an agency of the United States Government. Neither the United States Government nor the United States Department of Energy, nor Battelle, nor any of their employees, nor any jurisdiction or organization that has cooperated in the development of these materials, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness or any information, apparatus, product, software, or process disclosed, or represents that its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or any agency thereof, or Battelle Memorial Institute. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or any agency thereof.

PACIFIC NORTHWEST NATIONAL LABORATORY operated by BATTELLE for the UNITED STATES DEPARTMENT OF ENERGY under Contract DE-AC05-76RL01830
