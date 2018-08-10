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
 * \file cimmino.cpp
 * \brief Construction of the augmented systems and initialization of MUMPS
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include <abcd.h>
#include <fstream>

/*!
 *  \brief Construct the augmented systems and initialize the MUMPS solver
 *
 *  Construct the augmented systems depending on the partitions (grouped in block diagonal
 *  matrix), distribute the slaves to masters and place them on nodes, finally initialize
 *  MUMPS solver.
 *
 */
void abcd::initializeDirectSolver()
{
    int *sym_perm;

    mpi::broadcast(comm, icntl[Controls::nbparts], 0);

    if(comm.size() > parallel_cg) {
        if(instance_type == 0) {
            if(inter_comm.rank() == 0 && instance_type == 0)
                LINFO << "Initializing MUMPS";
            initializeMumps(mumps, true);
            createAugmentedSystems(n_aug, nz_aug, irn_aug, jcn_aug, val_aug);

            // initialize matrix in MUMPS
            mumps.n = n_aug;
            mumps.nz = nz_aug;
            mumps.irn = &irn_aug[0];
            mumps.jcn = &jcn_aug[0];
            mumps.a = &val_aug[0];

            if(inter_comm.rank() == 0 && instance_type == 0)
                LINFO << "Launching Initial MUMPS analysis";
            analyseAugmentedSystems(mumps);

            sym_perm = new int[mumps.n];
            std::copy(mumps.sym_perm, mumps.sym_perm + mumps.n, sym_perm);
        }

        allocateMumpsSlaves(mumps);

	// If the process is in a group with slaves (as master or not)
	if (my_slaves.size() > 0 || instance_type != 0) {
		// Then
		// As a master, we have to reinitialize mumps
		// As a slave, just initialize
	        if(instance_type == 0) {
	            mumps.job = -2;
	            dmumps_c(&mumps);
	            mumps.initialized = false;
	        }

	        initializeMumps(mumps);

                // Save ordering for masters with slaves
	        if(instance_type == 0) {
		    mumps.perm_in = new int[n_aug];
		    std::copy(sym_perm, sym_perm + n_aug, mumps.perm_in);
		    mumps.setIcntl(7, 1);
		}
	} else {
		// Useless participation in the blocking sub-communicator creation in
		// initializeMumps launched for Masters with slave..
		mpi::group grp;
		mpi::communicator(comm, grp);
	}
    } else {
        allocateMumpsSlaves(mumps);
        initializeMumps(mumps);
        createAugmentedSystems(n_aug, nz_aug, irn_aug, jcn_aug, val_aug);
    }

    if(instance_type == 0) {
        mumps.n = n_aug;
        mumps.nz = nz_aug;
        mumps.irn = &irn_aug[0];
        mumps.jcn = &jcn_aug[0];
        mumps.a = &val_aug[0];
    }
}    /* ----- end of method abcd::initializeDirectSolver ----- */
