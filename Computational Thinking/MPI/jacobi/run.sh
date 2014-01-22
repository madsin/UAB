#!/bin/bash

if [ $# -ne 2 ]
then
    echo "./run.sh numtasks c"
    exit -1
fi

if [ ! -f jacobi ]
then
    echo "Executable not found"
    exit -1
fi

openmpirun -np $1 jacobi $2
