########################
### Configuring qsub ###
########################

#$ -N MPI_measures
#$ -pe mpi-RR 16 			// Parallel environment; mpi-RR for round robin; mpi-FU for fill up
					// RR: divide processes among as many nodes as possible; FU: fill up nodes
#$ -q test.q				// Which queue to use; test.q for testing; aolin.q for actual measurements
#$ -v SGE_QMASTER_PORT
#$ -cwd
#$ -l h_rt=3600			        // Time limit

source /etc/profile.d/modules.sh	// Use this if module command cant be found

module load openmpi/1.6.3
