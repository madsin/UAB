#!/bin/bash

if [ $# -ne 4 ]
then
    echo "./run.sh numtasks dimM dimN dimO"
    exit -1
fi

if [ ! -f matmul_fancy ]
then
    echo "Executable not found"
    exit -1
fi

mpirun -np $1 matmul_fancy $2 $3 $4
