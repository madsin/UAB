#include <mpi.h>
#include <iostream>
#include <string>

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
    std::cout << "Hello from Process " << rank << " of " << numTasks << "!" << std::endl;

    /* Minimum 2 Processes */
    if ( numTasks < 2 ) {
        std::cout << "Minimum of 2 processes required" << std::endl;
        MPI_Abort ( MPI_COMM_WORLD, -1 );
        return -1;
    }

    if ( rank == 0 ) {
        std::string msg = "Hello from process 0!";
        int tag = 0;

        for ( int dest = 1; dest < numTasks; dest++ ) {
            retVal = MPI_Send ( (void *) buf, msg.size(), MPI_CHAR, dest, tag, MPI_COMM_WORLD );
        }
    }

    if ( rank != 0 ) {
    	char buf[256];
    	// TODO
    }

    /* Finalize MPI task */
    retVal = MPI_Finalize();

	return 0;
}
