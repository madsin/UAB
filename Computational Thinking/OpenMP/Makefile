GCC=gcc-mp-4.8 -Ofast -g -fopenmp -std=c99

default: all

all: original v3

icc_parallel:

original:
	export OMP_NUM_THREADS=1
	$(GCC) masterOMP_original.c -o gcc_masterOMP_original.out

v1:
	export OMP_NUM_THREADS=1
	$(GCC) masterOMP_v1.c -o gcc_masterOMP_v1.out

v2:
	export OMP_NUM_THREADS=8
	$(GCC) masterOMP_v2.c -o gcc_masterOMP_v2.out

v3:
	export OMP_NUM_THREADS=8
	$(GCC) masterOMP_v3.c -o gcc_masterOMP_v3.out

clean:
	rm -rf *.out *.dSYM

