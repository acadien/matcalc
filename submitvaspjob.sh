#!/bin/csh
#PBS -j oe
#PBS -N test
#PBS -l select=1:ncpus=8:mpiprocs=1
#PBS -l walltime=48:00:00

limit stacksize unlimited

setenv MACHTYPE x86-suse-linux
setenv TERM xterm
setenv MPI_BUFFER_MAX 2000000
setenv MPI_BUFS_PER_PROC 1024
setenv DAPL_MAX_CM_RESPONSE_TIME 22
setenv I_MPI_DEVICE rdssm
setenv I_MPI_DAPL_PROVIDER OpenIB-cma
setenv I_MPI_PIN_PROCS allcores
setenv I_MPI_PIN_MODE lib
setenv I_MPI_DEVICE_FALLBACK 0

echo "The nodefile is:"
cat $PBS_NODEFILE

setenv PROCS `cat $PBS_NODEFILE | wc -l`
setenv NODES `cat $PBS_NODEFILE | sort -u | wc -l`
echo "Running on $PROCS machines with 8 cores each."

cd ~asfd
echo "Working in directory: `pwd`"

mpdboot -r ssh -f $PBS_NODEFILE -n $NODES
mpiexec -np 12 ~echo > output 
mpdallexit

cd -
whoami >> /home/acadien/mozart_backup
echo ~asfd >> /home/acadien/mozart_backup
