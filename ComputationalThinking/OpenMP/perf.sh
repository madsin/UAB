#!/bin/bash

if [ $# -eq 0 ]
then
	echo "Need problem size as argument"
	exit -1
fi

# Measure performance of all existent executables
for EXEC in $(ls *out)
do

echo "$EXEC:"
perf stat -o ./$EXEC.log ./$EXEC $1

done
