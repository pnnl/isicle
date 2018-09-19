#!/bin/bash

module purge
export LD_LIBRARY_PATH=/usr/local/lib
module load intel/14.0.2
module load intelmpi

export ARMCI_DEFAULT_SHMMAX=32768
export NWCHEM_BASIS_LIBRARY="/people/scicons/apps/nwchem-6.6//src/basis/libraries/"
export NWCHEM_NWPW_LIBRARY="/people/scicons/apps/nwchem-6.6//src/nwpw/libraryps/"
# this disables xeon phi offload
export NWC_RANKS_PER_DEVICE=0
# this disables threaded in MKL since it is better to keep it to advanced users
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NWC_RANKS_PER_DEVICE=0
export ARMCI_OPENIB_DEVICE=mlx4_0
export OFFLOAD_INIT=on_offload
export MPIRETURN=999

FILE=$1

srun --mpi=pmi2 /people/scicons/apps/nwchem-6.6//bin/LINUX64/nwchem "$FILE" > "${FILE%.*}".out
