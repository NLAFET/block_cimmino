/*!
 * \file spldlt/init.cpp
 * \brief Initialization of SpLDLT solver
 * \author S. Cayrols
 * \version 1.0
 */

#include<abcd.h>
#include <spldlt.h>

/*!
 *  \brief Initialize SpLDLT parameters and structure
 *
 *  Initialize SpLDLT parameters and structure with the choice for verbosity
 *  and some features turned on.
 *
 *  \param spldlt: SPLDLT object
 *  \param local: first initialization for masters is local
 *
 */
void abcd::initializeSpLDLT(SPLDLT &spldlt, bool local)
{
  /* Initialize SpLDLT parameters and structure */

  spldlt.akeep = NULL;
  spldlt.fkeep = NULL;
  spldlt.ncpu  = 1;
  spldlt.ngpu  = 0;

  spldlt_default_options(&spldlt.options);

  spldlt.options.options.ordering    = 1;//use Metis ordering                          
  spldlt.options.options.scaling     = 0;//no scaling                                  
  spldlt.options.options.print_level = 1;//enable printing                             
  spldlt.options.options.print_level = 0;//disable printing                            
  spldlt.options.options.use_gpu     = (spldlt.ngpu > 0);//disable GPU

  spldlt_init(spldlt.ncpu, spldlt.ngpu);

}               /* -----  end of function abcd::initializeSpldlt  ----- */
