// Copyright Institut National Polytechnique de Toulouse (2014) 
// Contributor(s) :
// M. Zenadi <mzenadi@enseeiht.fr>
// D. Ruiz <ruiz@enseeiht.fr>
// R. Guivarch <guivarch@enseeiht.fr>

// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html"

// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 

// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 

// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.

/*!
 * \file mumps/analysis.cpp
 * \brief Analysis of MUMPS solver
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include <abcd.h>
#include <mumps.h>

/*!
 *  \brief Launch MUMPS analysis and display some info
 *
 *  Launch MUMPS analysis and display some info
 *
 *  \param mu: MUMPS object
 *
 */
void abcd::analyseAugmentedSystems(MUMPS &mu)
{

  double t = MPI_Wtime();

  // Run MUMPS Analysis
  mu(1);

  t = MPI_Wtime() - t;

  // Check if MUMPS succeded
  if(mu.getInfo(1) != 0){
    LERROR << string(32, '-') ;
    LERROR << "| MUMPS ANALYSIS FAILED on MA " << setw(7) << inter_comm.rank() << " |" ;
    LERROR << string(32, '-') ;
    LERROR << "| info(1)       : " << setw(6) << mu.getInfo(1) << string(4, ' ') << " |" ;
    LERROR << "| info(2)       : " << setw(6) << mu.getInfo(2) << string(4, ' ') << " |" ;
    LERROR << string(32, '-') ;;

    LERROR << "MUMPS exited with " << mumps_S.info[0];
    int job = -70 + mu.info[0];
    mpi::broadcast(intra_comm, job, 0);
    throw std::runtime_error("MUMPS exited with an error");
  }

  // Info display
  if(instance_type == 0) {
    double flop = mu.getRinfo(1);
    int prec = cout.precision();
    cout.precision(2);
    LINFO << string(32, '-') ;
    LINFO << "| MUMPS ANALYSIS on MA " << setw(7) << inter_comm.rank() << " |" ;
    LINFO << string(32, '-') ;
    LINFO << "| Flops estimate: " << setw(6) << scientific << flop << string(4, ' ') << " |" ;
    LINFO << "| Time:           " << setw(6) << t << " sec |" ;
    LINFO << string(32, '-') ;;
    cout.precision(prec);
  }
}               /* -----  end of function abcd::analyseAugmentedSystems  ----- */
