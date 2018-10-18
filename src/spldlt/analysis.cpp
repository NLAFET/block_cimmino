/*!
 * \file spldlt/analysis.cpp
 * \brief Analysis of SpLDLT solver
 * \author S. Cayrols
 * \version 1.0
 */

#include<abcd.h>
#include <spldlt.h>

/*!
 *  \brief Launch SpLDLT analysis and display some info
 *
 *  Launch SpLDLT analysis and display some info
 *
 *  \param inner_solver: SpLDLT object
 *
 */
void abcd::analyseAugmentedSystems(SPLDLT &inner_solver)
{

  double t = MPI_Wtime();

  // Run SpLDLT Analysis
  //TODO add the analysis here
  // spldlt_analyse();

}   /* -----  end of function abcd::analyseAugmentedSystems  ----- */
