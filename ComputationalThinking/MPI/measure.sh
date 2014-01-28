#!/bin/bash

# 10 minutes and test.q for now TODO
export QUEUE="aolin.q"
export QSUB0="qsub -q ${QUEUE} -v SGE_QMASTER_PORT -cwd -l h_rt=900"

export MPI_DIR=/home/master/master29/UAB/ComputationalThinking/MPI

# Cleanup previous logs
cd $MPI_DIR

# Matrix multiplication
echo "#########################"
echo "# Matrix Multiplication #"
echo "#########################"

for DIM in 2048 4096
do
    for N_TASKS in 8 16 32
    do
        cd $MPI_DIR/matmul
      	QSUB="$QSUB0 -N matmul_${N_TASKS}_${DIM}.log -pe mpi-RR $N_TASKS run.sh $N_TASKS $DIM $DIM"
        echo $QSUB
        echo $($QSUB)
        echo "#########################"
    done
done

# Matrix multiplication
echo "#########################"
echo "#        Jacobi         #"
echo "#########################"

# Jacobi
for C in 2 3
do
    for N_TASKS in 8 16 32
    do
        cd $MPI_DIR/jacobi
        QSUB="$QSUB0 -N jacobi_${N_TASKS}_${C}.log -pe mpi-RR $N_TASKS run.sh $N_TASKS $C"
        echo $QSUB
        echo $($QSUB)
        echo "#########################"
    done
done

exit 0
