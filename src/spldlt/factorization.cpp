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
  double t = MPI_Wtime();
  int posdef = 0;

  // Run SpLDLT factorization
  spldlt_factor(posdef, inner_solver.val, inner_solver.akeep,
      &inner_solver.fkeep, &inner_solver.options, &inner_solver.info);

}   /* -----  end of function abcd::factorizeAugmentedSystems  ----- */
