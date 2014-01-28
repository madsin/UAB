#!/bin/bash

if [ $# -ne 3 ]
then
    echo "./run.sh numtasks dimM dimN"
    exit -1
fi

if [ ! -f matmul ]
then
    echo "Executable not found"
    exit -1
fi

echo $(mpirun -np $1 matmul $2 $3)
