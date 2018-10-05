#!/bin/bash

source /etc/bashrc
module purge
module load nwchem/6.8.1_rhel7

cd /scratch

export ARMCI_DEFAULT_SHMMAX=131072
export NWCHEM_BASIS_LIBRARY="/home/scicons/cascade/apps/nwchem-6.8.1_rhel7/src/basis/libraries/"
export NWCHEM_NWPW_LIBRARY="/home/scicons/cascade/apps/nwchem-6.8.1_rhel7/src/nwpw/libraryps/"
#this disables xeon phi offload
export NWC_RANKS_PER_DEVICE=0
#this disables threaded in MKL since it is better to keep it to advanced users
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NWC_RANKS_PER_DEVICE=0
export ARMCI_OPENIB_DEVICE=mlx4_0
export OFFLOAD_INIT=on_offload

FILE=$1

srun --mpi=pmi2 -n ${SLURM_CPUS_ON_NODE} /dtemp/scicons/bin/nwchem6.8.1_rhel7 "$FILE" > "${FILE%.*}".out
