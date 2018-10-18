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

    // Run SpLDLT factorization
    //TODO add factorization of the augmented block
    // spldlt_factor()

}   /* -----  end of function abcd::factorizeAugmentedSystems  ----- */
