#!/bin/bash
#MSUB -l nodes=__:ppn=__,walltime=__,resfailpolicy=ignore
#MSUB -A sAccount
#MSUB -o __.%j.out
#MSUB -e __.%j.err
#MSUB -N g2m_00001
#MSUB -V
. /etc/profile.d/modules.sh
module purge
. /msc/apps/compilers/intel/15.0.090/composer_xe_2015.0.090/bin/compilervars.sh intel64
. /msc/apps/compilers/intel/impi/5.0.1.035/intel64/bin/mpivars.sh intel64

export ARMCI_DEFAULT_SHMMAX=32768
export NWCHEM_BASIS_LIBRARY="/home/scicons/cascade/apps/nwchem-6.6//src/basis/libraries/"
export NWCHEM_NWPW_LIBRARY="/home/scicons/cascade/apps/nwchem-6.6//src/nwpw/libraryps/"
#this disables xeon phi offload
export NWC_RANKS_PER_DEVICE=0
#this disables threaded in MKL since it is better to keep it to advanced users
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NWC_RANKS_PER_DEVICE=0
export ARMCI_OPENIB_DEVICE=mlx4_0
export OFFLOAD_INIT=on_offload
export MPIRETURN=999

cd __
MYDIR=$(/bin/pwd)
MYLIST=$(ls *.nw)
cd /dev/shm
for x in ${MYLIST}
do
 filename=`echo $x | sed -e 's@.nw@@'g`
 if [ ! -e $MYDIR/"$filename".output ]; then echo $MYDIR/"$filename".output
 srun -N 1 --exclusive --mpi=pmi2 -n __ /dtemp/scicons/bin/nwchem6.6 $MYDIR/"$x" > $MYDIR/"$filename".output &
 sleep 1
 fi
done
wait
export MPIRETURN=$?

############################################################################
# End of the job script
############################################################################

exit $MPIRETURN
