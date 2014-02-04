#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>

#define DEBUG 0
#define ROWS_PER_THREAD_C  (rowsPerThreadA)
#define COLS_PER_THREAD_C  (rowsPerThreadBt)

/******************************************************************************/
/* Matrix Multiplication: A (dimM x dimN) * B (dimN x dimO) = C (dimM x dimO) */
/******************************************************************************/

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
    int numTasks, rank;
    int dimM, dimN, dimO;
    int rowsPerThreadA, rowsPerThreadBt;
    int threadsPerRowCol;
    double *A, *Bt, *C, *myA, *myBt, *myC;
    double times[4];

    /* Initialize MPI task */
    retVal = MPI_Init ( &argc, &argv );
    if ( retVal != MPI_SUCCESS ) {
        std::cout << "Error initializing MPI" << std::endl;
        return -1;
    }

    MPI_Comm_size ( MPI_COMM_WORLD, &numTasks );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    /* Requests and Stats */
    MPI_Request reqsA[numTasks], reqsBt[numTasks];
    MPI_Status  statsA[numTasks], statsBt[numTasks];

    /* Parse arguments */
    if ( rank == 0 ) {
		if ( argc != 4 ) {
			std::cout << "Usage: ./matul dimM dimN dimO" << std::endl;
			MPI_Abort ( MPI_COMM_WORLD, -1 );
			return -1;
		} else {
			dimM = atoi(argv[1]);
			dimN = atoi(argv[2]);
			dimO = atoi(argv[3]);
			if (DEBUG) std::cout << "dimM = " << dimM
					             << ", dimN = " << dimN
					             << ", dimO = " << dimO << std::endl;
		}
    }
    /* Broadcast dimensions */
    MPI_Bcast ( &dimM, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast ( &dimN, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast ( &dimO, 1, MPI_INT, 0, MPI_COMM_WORLD );

    /***********************************/
    /* Check for valid data dimensions */
    /***********************************/

    /* Number of tasks must be a perfect square */
    if ( rank == 0 ) {
		if ( sqrt((double) numTasks) != floor(sqrt((double) numTasks)) ) {
			std::cout << "Number of tasks must be a perfect square" << std::endl;
			MPI_Abort ( MPI_COMM_WORLD, -1 );
			return -1;
		}
		threadsPerRowCol = (int) sqrt ((double) numTasks);
		if ( (dimM % threadsPerRowCol) || (dimO % threadsPerRowCol) ) {
			std::cout << "Dimensions not evenly dividable" << std::endl;
			MPI_Abort ( MPI_COMM_WORLD, -1 );
			return -1;
		}
    }
    MPI_Bcast ( &threadsPerRowCol, 1, MPI_INT, 0, MPI_COMM_WORLD );

    /******************/
    /* Initialization */
    /******************/

    /* Rows per thread */
    rowsPerThreadA  = dimM / threadsPerRowCol;
    rowsPerThreadBt = dimO / threadsPerRowCol;
    if ( (DEBUG) && (rank == 0) ) {
    	std::cout << "Each process will compute " << ROWS_PER_THREAD_C * COLS_PER_THREAD_C
    			  << " elements (gathered among " << rowsPerThreadA << " rows of A, "
    			  << rowsPerThreadBt << " rows of Bt (columns of B)" << std::endl;
    }

    myA  = (double *) malloc ( rowsPerThreadA * dimN * sizeof(double) );
    myBt = (double *) malloc ( rowsPerThreadBt * dimN * sizeof(double) );
    myC  = (double *) malloc ( ROWS_PER_THREAD_C * COLS_PER_THREAD_C * sizeof(double) );

    if ( rank == 0 ) {
        A  = (double *) malloc ( dimM * dimN * sizeof(double) );
        Bt = (double *) malloc ( dimO * dimN * sizeof(double) );
        C  = (double *) malloc ( dimM * dimO * sizeof(double) );

        srand ( time(NULL) );
        /* A, B random double in [0.5, 1.5], C = 0 */
        for ( int i = 0; i < dimM; i++ ) {
        	for ( int j = 0; j < dimN; j++ ) {
        		A [i*dimN + j] = 0.5 + ( (double) rand() / RAND_MAX );
        	}
        }
        for ( int i = 0; i < dimO; i++ ) {
        	for ( int j = 0; j < dimN; j++ ) {
        		Bt[i*dimN + j] = 0.5 + ( (double) rand() / RAND_MAX );
        	}
        }
        for ( int i = 0; i < dimM; i++ ) {
        	for ( int j = 0; j < dimO; j++ ) {
        		C[i*dimO + j] = 0;
        	}
        }

        if (DEBUG && dimM<12) {
			std::cout << "A=" << std::endl;
			printMat ( A, dimM, dimN );
			std::cout << "Bt=" << std::endl;
			printMat ( Bt, dimO, dimN );
    	}
    }
    for ( int i = 0; i < rowsPerThreadA*dimN; i++ ) {
    	myA[i] = 0;
    }
    for ( int i = 0; i < rowsPerThreadBt*dimN; i++ ) {
    	myBt[i] = 0;
    }
    for ( int i = 0; i < ROWS_PER_THREAD_C*COLS_PER_THREAD_C; i++ ) {
    	myC[i] = 0;
    }

    /* Take time */
    times[0] = MPI_Wtime();

    /***************/
    /* Divide data */
    /***************/

    /* Divide A and Bt */
    if ( rank == 0 ) {
    	int start;
    	if (DEBUG) std::cout << "P" << rank << ": MPI_ISend(), Time=" << MPI_Wtime()-times[0] << std::endl;
    	for ( int dest = 0; dest < numTasks; dest++ ) {
    		/* Calculate start address */
    		start = rowsPerThreadA*dimN*floor((double) dest/threadsPerRowCol);

    		/* Get own myA, myBt */
    		if ( dest == 0) {
				memcpy ( myA, &A[start], rowsPerThreadA*dimN*sizeof(double));
				memcpy ( myBt, &Bt[start], rowsPerThreadBt*dimN*sizeof(double));
				continue;
    		}

    		/* Non-blocking send */
    		MPI_Send ( &A[start], rowsPerThreadA*dimN, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD );
    		//MPI_Isend ( &A[start], rowsPerThreadA*dimN, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD, &reqsA[dest] );

    		start = rowsPerThreadBt*dimN*(dest % threadsPerRowCol);
    		/* Non-blocking send */
    		MPI_Send ( &Bt[start], rowsPerThreadBt*dimN, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD );
    		//MPI_Isend ( &Bt[start], rowsPerThreadBt*dimN, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD, &reqsBt[dest] );
    	}
    } else {
		/* Blocking receive */
		if (DEBUG) std::cout << "P" << rank << ": MPI_Recv(A), Time=" << MPI_Wtime()-times[0] << std::endl;
		MPI_Recv ( myA,  rowsPerThreadA*dimN,  MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &statsA[rank] );
		MPI_Recv ( myBt, rowsPerThreadBt*dimN, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &statsBt[rank] );
    }

	/* Debugging stuff */
    if (DEBUG && dimM<12) {
    	std::cout << "P" << rank << ": myA = ";
    	for ( int i = 0; i < rowsPerThreadA; i++ ) {
    		for ( int j = 0; j < dimN; j++ ) {
    			std::cout << myA[i*dimN + j] << " ";
    		}
    		std::cout << std::endl;
    	}
    }
    if (DEBUG && dimM<12) {
    	std::cout << "P" << rank << ": myBt = ";
    	for ( int i = 0; i < rowsPerThreadBt; i++ ) {
    		for ( int j = 0; j < dimN; j++ ) {
    			std::cout << myBt[i*dimN + j] << " ";
    		}
    		std::cout << std::endl;
    	}
    }

    /* Take time */
    times[1] = MPI_Wtime();

    /* Calculation */
    if (DEBUG) std::cout << "P" << rank << ": Calculation(), Time=" << MPI_Wtime()-times[0] << std::endl;
    for ( int i = 0; i < ROWS_PER_THREAD_C; i++ ) {
    	for ( int j = 0; j < COLS_PER_THREAD_C; j++ ) {
    		for ( int k = 0; k < dimN; k++ ) {
    			myC[i*COLS_PER_THREAD_C + j] += myA[i*dimN + k] * myBt[j*dimN + k];
    		}
    	}
    }

    /* Take time */
    times[2] = MPI_Wtime();

    if (DEBUG && dimM<12) {
    	std::cout << "P" << rank <<": myC = ";
    	for ( int i = 0; i < ROWS_PER_THREAD_C*COLS_PER_THREAD_C; i++ ) {
    		std::cout << myC[i] << " ";
    	}
    	std::cout << std::endl;
    }

	/* Data gathering */
    if (DEBUG) std::cout << "P" << rank << ": MPI_Gather(C), Time=" << MPI_Wtime()-times[0] << std::endl;
    MPI_Gather ( myC, ROWS_PER_THREAD_C*COLS_PER_THREAD_C, MPI_DOUBLE,
    		     C, ROWS_PER_THREAD_C*COLS_PER_THREAD_C, MPI_DOUBLE,
    		     0, MPI_COMM_WORLD );

    /* Take time */
    times[3] = MPI_Wtime();

    /* Data restructuring */
    double *Cr;
    if ( rank == 0) {
		Cr = (double *) malloc ( dimM * dimO * sizeof(double) );
		int m = 0;
		/* Iterate over submatrices */
		for ( int i = 0; i < threadsPerRowCol; i++ ) {
			for ( int j = 0; j < threadsPerRowCol; j++ ) {
				/* Iterate over elements of current submatrix */
				for ( int k = 0; k < ROWS_PER_THREAD_C; k++ ) {
					for ( int l = 0; l < COLS_PER_THREAD_C; l++ ) {
						Cr[dimO*(k + i*ROWS_PER_THREAD_C) + j*COLS_PER_THREAD_C + l] =C[m++];
					}
				}
			}
		}
    }

	/* Print results */
    if ( DEBUG && (rank == 0) && dimM<12 ) {
		std::cout << "C=" << std::endl;
		printMat ( C, dimM, dimO );
		std::cout << "Cr=" << std::endl;
		printMat ( Cr, dimM, dimO );
    }

	/* Print results */
    if ( rank == 0 ) std::cout << "RANK\tN_TASKS\tDIM\tEXEC\t\tCALC" << std::endl;
    MPI_Barrier ( MPI_COMM_WORLD );
    for ( int i = 0; i < numTasks; i++ ) {
    	if ( rank == i ) {
    		std::cout << rank << "\t" << numTasks << "\t" << dimN << "\t" << times[3]-times[0] << "\t\t" << times[2]-times[1] << std::endl;
    	}
    }

    /* Finalize MPI task */
    if (DEBUG) std::cout << "P" << rank << ": MPI_Finalize()" << std::endl;
    retVal = MPI_Finalize();

	/* Free memory */
    free(myA); free(myBt); free(myC);
	if ( rank == 0 ) {
		free(A); free(C); free(Bt);
	}

    return 0;
}

