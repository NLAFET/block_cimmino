#include "abcd.h"

int main(int argc, char* argv[])
{
  int err = 0;
  double t = 0.0;

  mpi::environment env(argc, argv);
  mpi::communicator world;
  string matrix_file, rhs_file, config_file;

  abcd solver;

  if(world.rank() == 0) {

    config_file = "config_file.info";
    if(argc >= 2) config_file = argv[1];

    std::cout << "Load the config_file " << config_file << std::endl;
    err = solver.parse_configFile(config_file, matrix_file, rhs_file);
    
    std::cout << "Load the matrix " << matrix_file << std::endl;
    err = solver.load_MM(matrix_file);

    if(!rhs_file.empty()){
      std::cout << "Load the RHS from " << rhs_file << std::endl;
      err = solver.load_RHS(rhs_file);
    }
  }

  try {
    t = MPI_Wtime();
    solver(BC_INITIALISE);
    solver(BC_PREPROCESS);
    solver(BC_FACTOR);
    solver(BC_SOLVE);
    if(world.rank() == 0) clog << "Total time: " << MPI_Wtime() - t << endl;
  } catch(std::runtime_error e) {
    std::cout << world.rank() << " Error code : " << e.what() << std::endl;
    err = 1;
  }
  world.barrier();

  return err;
}
