ISiCLE
======

Overview
--------
ISiCLE, or the _in silico_ chemical library engine, is a software framework for high-accuracy chemical property calculation. ISiCLE takes an [InChI](https://en.wikipedia.org/wiki/International_Chemical_Identifier) or [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) string as input, generates an initial 3D conformation, and subsequently optimizes this initial structure through molecular dynamics simulations and quantum chemistry optimizations. Finally, ISiCLE simulates desired properties (e.g. collision cross section, nuclear magnetic resonance chemical shifts) for each conformer yielded during molecular dynamics simulations to produce a single value, Boltzmann-weighted by relative Gibb's free energy, giving emphasis to properties from highly probable conformations.

<p align="center">
  <img align="center" src="docs/source/assets/schematic.svg" width="40%" height="40%">
</p>

Citing ISiCLE
-------------
If you would like to reference ISiCLE in an academic paper, we ask you include the following references:

* Colby, S.M., Thomas, D.G., Nuñez, J.R., Baxter, D.J., Glaesemann, K.R., Brown, J.M., Pirrung, M.A., Govind, N., Teeguarden, J.G., Metz, T.O. and Renslow, R.S., 2019. ISiCLE: A quantum chemistry pipeline for establishing in silico collision cross section libraries. _Analytical Chemistry_.
* Yesiltepe, Y., Nunez, J.R., Colby, S.M., Thomas, D.G., Borkum, M.I., Reardon, P.N., Washton, N.M., Metz, T.O., Teeguarden, J.T., Govind, N., and Renslow, R.S., 2018. An automated framework for NMR chemical shift calculations of small organic molecules. _Journal of Cheminformatics_.
* ISiCLE, version 2.0.0 http://github.com/pnnl/isicle (accessed MMM YYYY)

The first describes ISiCLE for CCS, the second describes ISiCLE for NMR chemical shifts, and the third is to cite the software package (update version and access date appropriately).

Disclaimer
----------
This material was prepared as an account of work sponsored by an agency of the United States Government.
Neither the United States Government nor the United States Department of Energy, nor Battelle, nor any of their employees, nor any jurisdiction or organization that has cooperated in the development of these materials, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness or any information, apparatus, product, software, or process disclosed, or represents that its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or any agency thereof, or Battelle Memorial Institute.
The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or any agency thereof.

PACIFIC NORTHWEST NATIONAL LABORATORY operated by BATTELLE for the UNITED STATES DEPARTMENT OF ENERGY under Contract DE-AC05-76RL01830
