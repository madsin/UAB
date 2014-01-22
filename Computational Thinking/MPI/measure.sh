#!/bin/bash

#$ -N MPI_measures
#$ -pe mpi 2
#$ -q test.q
#  -q aolin.q
#$ -v SGE_QMASTER_PORT
#$ -cwd
#$ -l h_rt=60

export MPI_DIR="/home/master/master29/UAB/Computational Thinking/MPI"
module load openmpi/1.6.3

# Matrix multiplication
for DIM in 512 1024 2048
do
    for N_TASKS in 2 4 8 16 32
    do
        CMD="bash $MPI_DIR/matmul/run.sh $N_TASKS $DIM $DIM"
        echo "Executing $CMD"
        #echo $($CMD)
    done
done

# Jacobi
for C in 1 2 4
do
    for N_TASKS in 2 4 8 16 32
    do
        CMD="bash $MPI_DIR/jacobi/run.sh $N_TASKS $C"
        echo "Executing $CMD"
        TEST="ls"
        #echo $($CMD)
    done
done

exit 0
