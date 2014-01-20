#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

#define SIZE 32

int main ( int argc, char **argv ) {
    int retVal;
    int numTasks, rank;

    /* Initialize MPI task */
    retVal = MPI_Init ( &argc, &argv );
    if ( retVal != MPI_SUCCESS ) {
        std::cout << "Error initializing MPI" << std::endl;
        return -1;
    }

    MPI_Comm_size ( MPI_COMM_WORLD, &numTasks );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    if ( rank == 0 ) std::cout << numTasks << "Processes" << std::endl;

    /* Minimum 2 Processes */
    if ( numTasks < 2 ) {
        std::cout << "Minimum of 2 processes required" << std::endl;
        MPI_Abort ( MPI_COMM_WORLD, -1 );
        return -1;
    }

    /* Initialize matrices */
    if ( rank == 0 ) {
        double *A, *Bt, *C;

        A  = (double *) malloc ( SIZE * SIZE * sizeof(double) );
        Bt = (double *) malloc ( SIZE * SIZE * sizeof(double) );
        C  = (double *) malloc ( SIZE * SIZE * sizeof(double) );

        srand ( time(NULL) );

        for ( int i = 0; i < SIZE; i++ ) {
        	for ( int j = 0; j < SIZE; j++ ) {
        		/* A, B random double in [0.5, 1.5] */
        		A [i*SIZE + j] = 0.5 + ( (double) rand() / RAND_MAX );
        		Bt[i*SIZE + j] = 0.5 + ( (double) rand() / RAND_MAX );
        		C [i*SIZE + j] = 0;
        	}
        }

        /* Scatter the data */
        /* Synchronous data divison */
        if ( !( SIZE % numTasks) ) {
        	int stepSize = SIZE / numTasks;
        	MPI_Request reqsA[stepSize * stepSize];
        	MPI_Request reqsB[stepSize * stepSize];
        	int itDest = 0;
        	int tag = 0;

        	/* Data can be divided into submatrices perfectly */
        	/* First submatrix is calculated by root process itself */
			for ( int itRow = stepSize; itRow < SIZE; itRow+=stepSize ) {
				for ( int itCol = stepSize; itCol < SIZE; itCol+=stepSize ) {
					/* Nonblocking send */
					itDest++;

					MPI_Isend ( (void *) &A[itRow * SIZE],
							    stepSize * SIZE,
							    MPI_DOUBLE,
							    itDest,
							    tag,
							    MPI_COMM_WORLD,
							    &reqsA[itDest] );

					MPI_Isend ( (void *) &Bt[itCol * SIZE],
								stepSize * SIZE,
								MPI_DOUBLE,
								itDest,
								tag,
								MPI_COMM_WORLD,
								&reqsB[itDest] );
				}
			}

		/* Asynchronous data divsion */
        } else {
        	/* TODO: Asynchronous data divison */
        	std::cout << "Asynchronous data divison not yet implemented" << std::endl;
        	return -1;
        }

        /* Do own calculation */

        /* Gather the data */

        /* Print results */

        /* Free memory */
        free(A); free(Bt); free(C);
    }

    /* Finalize MPI task */
    retVal = MPI_Finalize();

    return 0;
}
