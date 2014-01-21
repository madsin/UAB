#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

#define DEBUG 1
#define dimM 2
#define dimN 3

void printMat ( double *const M, int size ) {
	for ( int i = 0; i < size; i++ ) {
		for ( int j = 0; j < size; j++ ) {
			std::cout << M[i*size + j] << " ";
		}
		std::cout << std::endl;
	}
}

int main ( int argc, char **argv ) {
    int retVal;
    int numTasks, rank, stepSize;
    double *A, *Bt, *C, *myA, *myC;

    /* Initialize MPI task */
    retVal = MPI_Init ( &argc, &argv );
    if ( retVal != MPI_SUCCESS ) {
        std::cout << "Error initializing MPI" << std::endl;
        return -1;
    }

    MPI_Comm_size ( MPI_COMM_WORLD, &numTasks );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    /* Right now, only synchronous data division implemented */
    if ( ( dimM % numTasks ) && ( rank == 0 ) ) {
    	std::cout << "dimM % numTasks != 0, aborting" << std::endl;
    	MPI_Abort( MPI_COMM_WORLD, -1 );
    	return -1;
    }

    /* Initialization */
    Bt  = (double *) malloc ( dimM * dimM * sizeof(double) );
    myA = (double *) malloc ( dimM * sizeof(double) );
    myC = (double *) malloc ( dimM * sizeof(double) );

    if ( rank == 0 ) {
        A  = (double *) malloc ( dimM * dimM * sizeof(double) );
        C  = (double *) malloc ( dimM * dimM * sizeof(double) );

        srand ( time(NULL) );
        for ( int i = 0; i < dimM; i++ ) {
        	for ( int j = 0; j < dimM; j++ ) {
        		/* A, B random double in [0.5, 1.5] */
        		A [i*dimM + j] = 0.5 + ( (double) rand() / RAND_MAX );
        		Bt[i*dimM + j] = 0.5 + ( (double) rand() / RAND_MAX );
        		C [i*dimM + j] = 0;
        	}
        }
    }
    for ( int i = 0; i < dimM; i++ ) {
    	myA[i] = 0;
    	myC[i] = 0;
    }

    /* Divide data */
    /* Broadcast B */
    if (DEBUG) std::cout << "P" << rank << ": MPI_Bcast()" << std::endl;
    MPI_Bcast ( Bt, dimM*dimM, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    /* Scatter A */
    if (DEBUG) std::cout << "P" << rank << ": MPI_Scatter()" << std::endl;
    MPI_Scatter ( A, dimM, MPI_DOUBLE, myA, dimM, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if ( (DEBUG) && ( dimM < 8 ) ) {
    	std::cout << "P" << rank << ": myA = ";
    	for ( int i = 0; i < dimM; i++ ) {
    		std::cout << myA[i] << " ";
    	}
    	std::cout << std::endl;
    }

    /* Calculation */
    if (DEBUG) std::cout << "P" << rank << ": Calculation()" << std::endl;
    for ( int i = 0; i < dimM; i++ ) {
    	for ( int j = 0; j < dimM; j++ ) {
    		myC[i] += myA[j] * Bt[i*dimM + j];
    	}
    }
    if ( (DEBUG) && ( dimM < 8 ) ) {
    	std::cout << "P" << rank <<": myC = ";
    	for ( int i = 0; i < dimM; i++ ) {
    		std::cout << myC[i] << " ";
    	}
    	std::cout << std::endl;
    }

	/* Data gathering */
    if (DEBUG) std::cout << "P" << rank << ": MPI_Gather()" << std::endl;
    MPI_Gather( myC, dimM, MPI_DOUBLE, C, dimM, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	/* Print results */
    if ( DEBUG && (rank == 0) && (dimM <= 8) ) {
		std::cout << "A=" << std::endl;
		printMat ( A, dimM );
		std::cout << "Bt=" << std::endl;
		printMat ( Bt, dimM );
		std::cout << "C=" << std::endl;
		printMat ( C, dimM );
    }

    /* Finalize MPI task */
    if (DEBUG) std::cout << "P" << rank << ": MPI_Finalize()" << std::endl;
    retVal = MPI_Finalize();

	/* Free memory */
    free(myA); free(Bt); free(myC);
	if ( rank == 0 ) {
		free(A); free(C);
	}

	/* Checksum calc */
	if ( rank == 0 ) {
		double checksum_parallel = 0;
		for ( int i = 0; i < dimM; i++ ) {
			for ( int j = 0; j < dimM; j++ ) {
				checksum_parallel += C[i*dimM + j];
			}
		}
		std::cout << "Parallel checksum result = " << checksum_parallel << std::endl;
	}

    return 0;
}
