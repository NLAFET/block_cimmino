#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include "abcd_c.h"

int main(int argc, char **argv)
{
  int err = 0;
  int rank = 0;
  double t = 0.0;
  char *matrix_file = NULL;
  char *rhs_file    = NULL;
  char config_file[50] = "config_file.info";
  abcd_c *solver;

  /* Initialize MPI */
  MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  solver = new_solver();

  if(rank == 0) { // the master

    if(argc >= 2) strncpy(config_file, argv[1], strlen(argv[1]) + 1);

    printf("Load the config_file %s\n", config_file);
    err = parse_configFile(solver, config_file, &matrix_file, &rhs_file);
    
    printf("Load the matrix %s\n", matrix_file);
    err = load_MM(solver, matrix_file);

    if(rhs_file){
      printf("Load the RHS from %s\n", rhs_file);
      err = load_RHS(solver, rhs_file);
    }
  }

  t = MPI_Wtime();
  call_solver(solver, BC_INITIALISE);
  call_solver(solver, BC_PREPROCESS);
  call_solver(solver, BC_FACTOR);
  call_solver(solver, BC_SOLVE);
  if(rank == 0) printf("Total time: %f\n", MPI_Wtime() - t);

  MPI_Finalize();

  return 0;
}
