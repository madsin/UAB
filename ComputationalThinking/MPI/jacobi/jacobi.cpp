#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>

#define DEBUG 0
#define DIM (1024)
#define PRECISION (1.0e-6)

void printMat ( double *M, int size ) {
	for ( int i = 0; i < size; i++ ) {
		for ( int j = 0; j < size; j++ ) {
			std::cout << M[i*size + j] << " ";
		}
		std::cout << std::endl;
	}
}

void printVec ( double *vec, int size ) {
	for ( int i = 0; i < size; i++ ) {
		std::cout << vec[i] << " ";
	}
	std::cout << std::endl;
}

int main ( int argc, char **argv ) {
    int retVal;
    int numTasks, rank, stepSize;
    int itCount = 0;
    int dim, c;
    double times[4], calcTime = 0;
    double epsilon = 1;
    double *A, *b, *x, *x_old, *myA, *myB, *myX_new;

    /* Initialize MPI task */
    retVal = MPI_Init ( &argc, &argv );
    if ( retVal != MPI_SUCCESS ) {
        std::cout << "Error initializing MPI" << std::endl;
        return -1;
    }
    MPI_Comm_size ( MPI_COMM_WORLD, &numTasks );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    /* Get dimension */
    if ( rank == 0 ) {
    	if ( argc != 2 ) {
    		std::cout << "Usage: ./jacobi c" << std::endl;
    		MPI_Abort ( MPI_COMM_WORLD, -1 );
    		return -1;
    	} else {
    		c = atoi(argv[1]);
    		dim = c*DIM;
    		if (DEBUG) std::cout << "dim = " << dim << std::endl;
    	}
    }
    MPI_Bcast ( &dim, 1, MPI_INT, 0, MPI_COMM_WORLD );

    /* Only synchronous data division supported */
    if ( (dim % numTasks) && (rank == 0) ) {
    	std::cout << "dim % numTasks != 0. Aborting." << std::endl;
    	MPI_Abort ( MPI_COMM_WORLD, -1 );
    	return -1;
    }
    stepSize = dim / numTasks;
    if ( (DEBUG) && (rank == 0) ) std::cout << "stepSize = " << stepSize << std::endl;

    /* Initialize */
    if (DEBUG) std::cout << "P" << rank << ": Init()" << std::endl;
    if ( rank == 0 ) {
    	A = (double *) malloc ( dim * dim * sizeof(double) );

        /* Create a diagonal dominant matrix A */
    	srand ( time(NULL) );
        for ( int i = 0; i < dim; i++ ) {
        	A[i*dim + i] = 0;
        	for ( int j = 0; j < dim; j++ ) {
        		if ( i == j ) continue;
        		A[i*dim + j] = 0.5 + ( (double) rand() / RAND_MAX );
        		A[i*dim + i] += A[i*dim + j];
        	}
        	A[i*dim + i] += 0.5 + ( (double) rand() / RAND_MAX );
        }

        /* Check if A diagonal dominant */
        if (DEBUG) {
        	for ( int i = 0; i < dim; i++ ) {
        		int sum = 0;
        		for ( int j = 0; j < dim; j++ ) {
        			if ( i == j ) continue;
        			sum += A[i*dim + j];
        		}
        		if ( A[i*dim + i] <= sum ) {
        			std::cout << "A is not diagonal dominant. Aborting." << std::endl;
        			MPI_Abort ( MPI_COMM_WORLD, -1 );
        			return -1;
        		}
        	}
        std::cout << "Created diagonal dominant matrix A=" << std::endl;
        printMat( A, dim );
        }

        /* Create vector b */
        b     = (double *) malloc ( dim * sizeof(double) );
        x     = (double *) malloc ( dim * sizeof(double) );
        x_old = (double *) malloc ( dim * sizeof(double) );
        for ( int i = 0; i < dim; i++ ) {
        	b[i] = 0.5 + ( (double) rand() / RAND_MAX );
        	// Initial solution
        	x[i] = 0.5 + ( (double) rand() / RAND_MAX );
        }
        if (DEBUG) {
        	std::cout << "b =" << std::endl;
        	printVec ( b, dim );
        	std::cout << "x_0 =" << std::endl;
        	printVec ( x, dim );
        }
    } else {
    	x = (double *) malloc ( dim * sizeof(double) );
    }
    myA     = (double *) malloc ( stepSize * dim * sizeof(double) );
    myB     = (double *) malloc ( stepSize * sizeof(double) );
    myX_new = (double *) malloc ( stepSize * sizeof(double) );

    /* Take time */
    times[0] = MPI_Wtime();

    /* Synchronous data distribution */
    // A
    if (DEBUG) std::cout << "P" << rank << ": MPI_Scatter(A)" << std::endl;
    MPI_Scatter ( A, stepSize*dim, MPI_DOUBLE,
    		      myA, stepSize*dim, MPI_DOUBLE,
    		      0, MPI_COMM_WORLD );
    if (DEBUG) {
    	std::cout << "P" << rank << ": myA=" << std::endl;
    	for ( int i = 0; i < stepSize; i++ ) {
    		printVec ( &myA[i*dim], dim );
    	}
    }
    // b
    if (DEBUG) std::cout << "P" << rank << ": MPI_Scatter(b)" << std::endl;
    MPI_Scatter ( b, stepSize, MPI_DOUBLE,
    		      myB, stepSize, MPI_DOUBLE,
    		      0, MPI_COMM_WORLD );
    if (DEBUG) {
    	std::cout << "P" << rank << ": myB=" << std::endl;
    	for ( int i = 0; i < stepSize; i++ ) {
    		std::cout << myB[i] << " ";
    	}
    	std::cout << std::endl;
    }

    /* Iterate */
    if (DEBUG) std::cout << "P" << rank << ": Iterate()" << std::endl;
    while ( epsilon > PRECISION ) {
    	itCount++;
    	if (DEBUG) {
    		std::cout << "P" << rank << ": Starting iteration " << itCount << std::endl;
    	}

    	/* Store old x */
    	if ( rank == 0 ) {
    		for ( int i = 0; i < dim; i++ ) {
    			x_old[i] = x[i];
    		}
    	}

    	/* Incremental synchronous data distribution */
    	if (DEBUG) std::cout << "P" << rank << ": MPI_Bcast(x)" << std::endl;
    	MPI_Bcast ( x, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    	if ( (rank == 0) && (DEBUG) ) {
    		std::cout << "x_" << itCount-1 << "=" << std::endl;
    		printVec ( x, dim );
    	}

    	/* Take time */
    	times[1] = MPI_Wtime();

    	/* Calculate */
    	if (DEBUG) std::cout << "P" << rank << ": Calculate()" << std::endl;
    	for ( int i = 0; i < stepSize; i++ ) {
    		myX_new[i] = myB[i];
    		for ( int j = 0; j < dim; j++ ) {
    			if ( (rank*stepSize + i) == j ) continue;

    			myX_new[i] -= myA[i*dim + j] * x[j];
    		}
    		myX_new[i] /= myA[i*dim + i + rank*stepSize];
    		if (DEBUG) std::cout << "P" << rank << ": myX_new[" << i << "]="
    				             << myX_new[i] << std::endl;
    	}

    	/* Take time */
    	times[2] = MPI_Wtime();
    	calcTime += times[2] - times[1];

    	/* Gather data */
    	if (DEBUG) std::cout << "P" << rank << ": MPI_Gather(x)" << std::endl;
    	MPI_Gather( myX_new, stepSize, MPI_DOUBLE,
    			    x, stepSize, MPI_DOUBLE,
    			    0, MPI_COMM_WORLD );
    	if ( (DEBUG) && (rank == 0) ) {
    		std::cout << "x_" << itCount << " =" << std::endl;
    		printVec ( x, dim );
    	}

    	/* Calculate epsilon */
    	if ( rank == 0 ) {
    		if (DEBUG) std::cout << "P" << rank << ": epsilon()" << std::endl;
    		epsilon = fabs ( x[0] - x_old[0] );
    		for ( int i = 1; i < dim; i++ ) {
    			double temp = fabs ( x[i] - x_old[i] );
    			if ( epsilon < temp ) epsilon = temp;
    		}
    	}
    	if (DEBUG) {
    		if (rank == 0) std::cout << "epsilon_" << itCount << " = " << epsilon << std::endl;
    		std::cout << "P" << rank << ": MPI_Bcast(epsilon)" << std::endl;
    	}
    	MPI_Bcast ( &epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    }

    /* Take time */
    times[3] = MPI_Wtime();

    /* Finalize MPI task */
    retVal = MPI_Finalize();

    /* Print results */
    if ( rank == 0 ) std::cout << "Total of " << itCount << " iterations performed" << std::endl;
    std::cout << "P" << rank << ": Execution time = " << (times[1] - times[0]) << std::endl
    		                 << "Calculation time = " << calcTime << std::endl;

    /* Free Memory */
    if ( rank == 0) {
    	free(A); free(x); free(b);
    }
    free(myA); free(myB); free(myX_new);

    return 0;
}
