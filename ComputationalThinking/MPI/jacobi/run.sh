#!/bin/bash

source /etc/profile.d/modules.sh
module load openmpi/1.6.3

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

echo $(mpirun -np $1 jacobi $2)
