#!/bin/bash

if [ $# -ne 4 ]
then
  echo "Usage: `basename $0` <command to run> <directory to run from> <# Procs [1,8,16..]> <name of run>"
  exit 65
fi

COMMAND=$1
DIRECTORY=$2
if [ ${DIRECTORY:0:1} != "/" ] && [ ${DIRECTORY:0:1} != "~" ]
then
    echo "Directory must be an absolute path (start with '/' or '~')."
    exit 65
fi
NPROCS=$3
NAME=$4

SEL=1
NCPUS=$NPROCS
if [ $NCPUS -gt 1 ]
then
  SEL=$(( $NPROCS / 8 ))
  NCPUS=8
fi

SUBFILE=submitvaspjob.sh

cat > $SUBFILE << EOF
#!/bin/csh
#PBS -j oe
#PBS -N $NAME
#PBS -l select=$SEL:ncpus=$NCPUS:mpiprocs=1
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
cat \$PBS_NODEFILE

setenv PROCS \`cat \$PBS_NODEFILE | wc -l\`
setenv NODES \`cat \$PBS_NODEFILE | sort -u | wc -l\`
echo "Running on \$PROCS machines with $NCPUS cores each."

cd $DIRECTORY
echo "Working in directory: \`pwd\`"

mpdboot -r ssh -f \$PBS_NODEFILE -n \$NODES
mpiexec -np $NPROCS $COMMAND > output 
mpdallexit

cd -
whoami >> /home/acadien/mozart_backup
echo $DIRECTORY >> /home/acadien/mozart_backup
EOF

qsub $SUBFILE