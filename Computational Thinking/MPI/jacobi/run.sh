#!/bin/bash

if [ $# -ne 1 ]
then
    echo "./run.sh numtasks"
    exit -1
fi

if [ ! -f jacobi ]
then
    echo "Executable not found"
    exit -1
fi

openmpirun -np $1 jacobi
