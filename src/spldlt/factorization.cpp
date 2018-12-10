/*!
 * \file mumps/factorization.cpp
 * \brief Factorization using SpLDLT solver
 * \author S. Cayrols
 * \version 1.0
 */

#include<abcd.h>
#include <spldlt.h>

/*!
 *  \brief Launch SpLDLT factorization and display some info
 *
 *  Launch SpLDLT factorization and display some info
 *
 *  \param inner_solver: SpLDLT object
 *
 */
void abcd::factorizeAugmentedSystems(SPLDLT &inner_solver)
{
  double t      = MPI_Wtime();
  int posdef    = 0;
  long flops    = 0;
  long nfactor  = 0;
  long mem      = 0;
  int error     = 0;
  int gerror    = 0;

  flops   = inner_solver.info.num_flops;
  nfactor = inner_solver.info.num_factor;
//mem     = inner_solver.info.factor_mem; //TODO ?

  // Run SpLDLT factorization
  spldlt_factor(posdef, inner_solver.val, inner_solver.akeep,
      &inner_solver.fkeep, &inner_solver.options, &inner_solver.info);

  t = MPI_Wtime() - t;

  if(inner_solver.info.flag != 0){
    error = 1;
  }

  mpi::all_reduce(inter_comm, error, gerror, std::plus<int>());
  //mpi::reduce(inter_comm, vol, totVol, std::plus<int>(),0);

  if(gerror > 0){
    LINFO << "SpLDLT returned error " << inner_solver.info.flag;
    throw std::runtime_error("Factorisation exited with an error!");
  }

  // Check if SpLDLT succeded
  if(inner_solver.info.flag != 0){
    LINFO << "SpLDLT returned error " << inner_solver.info.flag;
    throw std::runtime_error("SpLDLT exited with an error!");
  }

  // Info display
  if(instance_type == 0) {

    int prec = cout.precision();
    cout.precision(2);
    LINFO << string(32, '-') ;
    LINFO << "| SpLDLT FACTORIZ on MA " << setw(7) << inter_comm.rank() << " |" ;
    LINFO << string(32, '-') ;
    LINFO << "| N             : " << setw(12) << inner_solver.n << " |" ;
    LINFO << "| NZ            : " << setw(12) << inner_solver.nz << " |" ;
    LINFO << "| FNZ           : " << setw(12) << nfactor << " |" ;
    LINFO << "| Flops         : " << setw(6) << scientific << flops << string(4, ' ') << " |" ;
    LINFO << "| Time          : " << setw(6) << t << " sec |" ;
    LINFO << "| Average memory    : " << setw(6) << mem << " M| ";
    LINFO << "| Delay             : " << setw(6) << inner_solver.info.num_delay << " |" ;
    LINFO << string(32, '-') ;;
    cout << "SpLDLT STATS " << inter_comm.rank() << " : " << inner_solver.n
      << " " << nfactor << endl;
    cout.precision(prec);
  }

  if(IRANK == 0) LINFO << "Factorization average memory : " << setw(6) 
    << nfactor << " M";
  if(IRANK == 0) LINFO << "Factorization maximum memory : " << setw(6) 
    << mem << " M";
  if(IRANK == 0) LINFO << "Factorization total flops : " << setw(6) 
    << flops << " flops";
}   /* -----  end of function abcd::factorizeAugmentedSystems  ----- */
