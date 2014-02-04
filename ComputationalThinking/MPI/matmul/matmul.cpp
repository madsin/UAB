#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

#define DEBUG 1

void printMat ( double *const M, const int sizeM, const int sizeN ) {
	for ( int i = 0; i < sizeM; i++ ) {
		for ( int j = 0; j < sizeN; j++ ) {
			std::cout << M[i*sizeN + j] << " ";
		}
		std::cout << std::endl;
	}
}

int main ( int argc, char **argv ) {
    int retVal;
    int numTasks, rank, stepSize;
    int dimM, dimN;
    double *A, *Bt, *C, *myA, *myC;
    double times[4];

    /* Initialize MPI task */
    retVal = MPI_Init ( &argc, &argv );
    if ( retVal != MPI_SUCCESS ) {
        std::cout << "Error initializing MPI" << std::endl;
        return -1;
    }

    MPI_Comm_size ( MPI_COMM_WORLD, &numTasks );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    /* Parse arguments */
    if ( rank == 0 ) {
		if ( argc != 3 ) {
			std::cout << "Usage: ./matul dimM dimN" << std::endl;
			MPI_Abort ( MPI_COMM_WORLD, -1 );
			return -1;
		} else {
			dimM = atoi(argv[1]);
			dimN = atoi(argv[2]);
			if (DEBUG) std::cout << "dimM = " << dimM
					           << ", dimN = " << dimN << std::endl;
		}
    }
    MPI_Bcast ( &dimM, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast ( &dimN, 1, MPI_INT, 0, MPI_COMM_WORLD );

    /* Right now, only synchronous data division implemented */
    if ( ( dimM % numTasks ) && ( rank == 0 ) ) {
    	std::cout << "dimM % numTasks != 0, aborting" << std::endl;
    	MPI_Abort( MPI_COMM_WORLD, -1 );
    	return -1;
    }

    /* Initialization */
    stepSize = dimM / numTasks;
    if ( (DEBUG) && (rank == 0) ) std::cout << "stepSize = " << stepSize << std::endl;

    Bt  = (double *) malloc ( dimM * dimN * sizeof(double) );
    myA = (double *) malloc ( stepSize * dimN * sizeof(double) );
    myC = (double *) malloc ( stepSize * dimM * sizeof(double) );

    if ( rank == 0 ) {
        A  = (double *) malloc ( dimM * dimN * sizeof(double) );
        C  = (double *) malloc ( dimM * dimM * sizeof(double) );

        srand ( time(NULL) );
        for ( int i = 0; i < dimM; i++ ) {
        	for ( int j = 0; j < dimN; j++ ) {
        		/* A, B random double in [0.5, 1.5] */
        		A [i*dimN + j] = 0.5 + ( (double) rand() / RAND_MAX );
        		Bt[j*dimM + i] = 0.5 + ( (double) rand() / RAND_MAX );
        	}
        	for ( int j = 0; j < dimM; j++ ) {
        		C[i*dimM + j] = 0;
        	}
        }

        if (DEBUG && dimM<12 ) {
			std::cout << "A=" << std::endl;
			printMat ( A, dimM, dimN );
			std::cout << "Bt=" << std::endl;
			printMat ( Bt, dimM, dimN );
    	}
    }
    for ( int i = 0; i < stepSize*dimN; i++ ) {
    	myA[i] = 0;
    }
    for ( int i = 0; i < stepSize*dimM; i++ ) {
    	myC[i] = 0;
    }

    /* Take time */
    times[0] = MPI_Wtime();

    /* Divide data */
    /* Broadcast B */
    MPI_Status bla[numTasks];
    if (rank == 0) {
    	for ( int i = 1; i < numTasks; i++ ) {
    		if (DEBUG) std::cout << "P" << rank << ": MPI_Send(Bt), Time=" << MPI_Wtime()-times[0] << std::endl;
    		MPI_Send ( Bt, dimM*dimN, MPI_DOUBLE, i, 0, MPI_COMM_WORLD );
    	}
    } else {
    	if (DEBUG) std::cout << "P" << rank << ": MPI_Recv(Bt), Time=" << MPI_Wtime()-times[0] << std::endl;
    	MPI_Recv( Bt, dimM*dimN, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &bla[rank] );
    }

    //if (DEBUG) std::cout << "P" << rank << ": MPI_Bcast(), Time=" << MPI_Wtime()-times[0] << std::endl;
    //MPI_Bcast ( Bt, dimM*dimN, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    /* Scatter A */
    if (DEBUG) std::cout << "P" << rank << ": MPI_Scatter(), Time=" << MPI_Wtime()-times[0] << std::endl;
    MPI_Scatter ( A, stepSize*dimN, MPI_DOUBLE, myA, stepSize*dimN, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if ( (DEBUG) && ( dimM < 12 ) ) {
    	std::cout << "P" << rank << ": myA = ";
    	for ( int i = 0; i < stepSize*dimN; i++ ) {
    		std::cout << myA[i] << " ";
    	}
    	std::cout << std::endl;
    }

    /* Take time */
    times[1] = MPI_Wtime();

    /* Calculation */
    if (DEBUG) std::cout << "P" << rank << ": Calculation(), Time=" << MPI_Wtime()-times[0] << std::endl;
    for ( int i = 0; i < stepSize; i++ ) {
    	for ( int j = 0; j < dimM; j++ ) {
    		for ( int k = 0; k < dimN; k++ ) {
    			myC[i*dimM + j] += myA[i*dimN + k] * Bt[(i+j)*dimN + k];
    		}
    	}
    }
    /* Take time */
    times[2] = MPI_Wtime();
    if ( (DEBUG) && ( dimM < 8 ) ) {
    	std::cout << "P" << rank <<": myC = ";
    	for ( int i = 0; i < stepSize*dimM; i++ ) {
    		std::cout << myC[i] << " ";
    	}
    	std::cout << std::endl;
    }

	/* Data gathering */
    if (DEBUG) std::cout << "P" << rank << ": MPI_Gather(), Time=" << MPI_Wtime()-times[0] << std::endl;
    MPI_Gather( myC, stepSize*dimM, MPI_DOUBLE, C, stepSize*dimM, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	/* Print results */
    if ( DEBUG && (rank == 0) && (dimM <= 8) ) {
		std::cout << "C=" << std::endl;
		printMat ( C, dimM, dimM );
    }

    /* Take time */
    times[3] = MPI_Wtime();

	/* Print results */
    if ( rank == 0 ) std::cout << "RANK\tN_TASKS\tDIM\tEXEC\t\tCALC" << std::endl;
    MPI_Barrier ( MPI_COMM_WORLD );
    for ( int i = 0; i < numTasks; i++ ) {
    	if ( rank == i ) {
    		std::cout << rank << "\t" << numTasks << "\t" << dimN << "\t" << times[3]-times[0] << "\t\t" << times[2]-times[1] << std::endl;
    	}
    }

    /* Finalize MPI task */
    if (DEBUG) std::cout << "P" << rank << ": MPI_Finalize(), Time=" << MPI_Wtime()-times[0] << std::endl;
    retVal = MPI_Finalize();

	/* Free memory */
    free(myA); free(Bt); free(myC);
	if ( rank == 0 ) {
		free(A); free(C);
	}

    return 0;
}
