#!/bin/bash

#$ -N MPI_measures
#$ -pe mpi 2
#$ -q test.q
#$ -v SGE_QMASTER_PORT
#$ -cwd
#$ -l h_rt=60
source /etc/profile.d/modules.sh
module load openmpi/1.6.3

export MPI_DIR=/home/master/master29/UAB/ComputationalThinking/MPI

# Matrix multiplication
for DIM in 512 #1024 2048
do
    for N_TASKS in 2 4 8 #16 32
    do
        echo "Executing /matmul/run.sh $N_TASKS $DIM $DIM"
        cd $MPI_DIR/matmul
        pwd
        bash ./run.sh $N_TASKS $DIM $DIM
    done
done

# Jacobi
for C in 1 #2 4
do
    for N_TASKS in 2 #4 8 16 32
    do
        echo "Executing /jacobi/run.sh $N_TASKS $C"
        cd $MPI_DIR/jacobi
        bash ./run.sh $N_TASKS $C
    done
done

exit 0
