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
 * \file worker.cpp
 * \brief Implementation of the waiting-to-help state of the slaves
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include <abcd.h>

/*!
 *  \brief wait to help the master with MUMPS
 *
 *  In an infinite loop:
 *  receive an order from the master (job) then wait until the process can help it
 *  with MUMPS. Orders can be:
 *   - -1: Go Home G.I.
 *   - 1: Solve augmented system (simple or not)
 *   - 2: solve S after initializing with arrays from job 3
 *   - 3: solve augmented systems with distributed output upon building of S
 *   - 4: solve S without re-building it
 *
 */
    void
abcd::waitForSolve()
{
    DMUMPS_STRUC_C mu;

    // arrays of rows/columns/values for S
    std::vector<int> vrows, vcols;
    std::vector<double> vvals;

    std::vector<int> target;
    std::vector<int> target_idx;
    target.reserve(size_c);
    target_idx.reserve(size_c);

    int job = 0;
    do{
        // receive order
        mpi::broadcast(intra_comm, job, 0);

        // if we receive an emergency stop
        if(job < -1){
            info[Controls::status] = job;
            throw runtime_error("Worker was asked to stop by its master");
        }

        //we are finished here!
        if(job == -1) break;
        // do the solve, in the (simple or not) sumproject
        if(job == 1){
            mumps(3);
        // handle the matrix S
        } else if (job == 2) {
            // initialize MUMPS for S
            mumps_S.sym = 2;
            mumps_S.par = 1;
            mumps_S.comm_fortran = MPI_Comm_c2f((MPI_Comm) comm);
            mumps_S(-1);

            // If MUMPS verbose chosen
            if (icntl[Controls::mumps_verbose]) {
                mumps_S.setIcntl(1, 6);
                mumps_S.setIcntl(2, 0);
                mumps_S.setIcntl(3, 6);
                mumps_S.setIcntl(4, 2);
            // or silent
            } else {
                mumps_S.setIcntl(1, -1);
                mumps_S.setIcntl(2, -1);
                mumps_S.setIcntl(3, -1);
                mumps_S.setIcntl(4, -1);
            }

            // use S from job 3 to initialize matrix in MUMPS
            mumps_S.nz_loc = vrows.size();
            mumps_S.irn_loc = &vrows[0];
            mumps_S.jcn_loc = &vcols[0];
            mumps_S.a_loc = &vvals[0];

            // Launch MUMPS analysis/factorization/solve on S
            mumps_S(1);
            if(mumps_S.info[0] < 0) continue;
            mumps_S(2);
            if(mumps_S.info[0] < 0) continue;
            mumps_S(3);
        // do the solve with a distributed solution output, when building S
        } else if (job == 3) {
            int s;
            int start_c;
            // get some necessary arrays and data
            mpi::broadcast(intra_comm, s, 0);
            std::vector<int> loc_cols, mycols;
            mpi::broadcast(intra_comm, start_c, 0);

            if(start_c == -1)
                continue;

            mpi::broadcast(intra_comm, column_index, 0);
            mpi::broadcast(intra_comm, loc_cols, 0);
            mpi::broadcast(intra_comm, n_o, 0);
            mpi::broadcast(intra_comm, mycols, 0);

            // initialize local solution
            mumps.lsol_loc = mumps.getInfo(23);
            mumps.sol_loc = new double[mumps.lsol_loc * s];
            mumps.isol_loc = new int[mumps.lsol_loc];

            mumps(3);
            if(mumps.info[0] < 0) continue;

            int part = 0;
            int end_c = column_index[part].size();
            int i_loc = 0;

            // compress local solution
            if(target.size() == 0)
                while (i_loc < mumps.lsol_loc){
                    int isol = mumps.isol_loc[i_loc] - 1;
                    if (isol >= start_c && isol < end_c){
                        int ci = column_index[part][isol] - n_o;
                        target.push_back(i_loc);
                        target_idx.push_back(ci);
                    }
                    i_loc++;
                }

            // initialize local part of the S matrix (0.5-Ai+Aiel)
            for (int j = 0; j < s; j++) {
                int col = mycols[j];
                for (int i = 0; i < (int)target_idx.size(); i++) {
                    if(target_idx[i] < col)
                        continue;
                    else if ( target_idx[i] == col)
                        vvals.push_back(0.5 - mumps.sol_loc[target[i] + j * mumps.lsol_loc]);
                    else
                        vvals.push_back(-mumps.sol_loc[target[i] + j * mumps.lsol_loc]);
                    vrows.push_back(target_idx[i] + 1);
                    vcols.push_back(col + 1);

                }
            }

            delete[] mumps.isol_loc;
            delete[] mumps.sol_loc;
        // handle the matrix S w/o re-building it
        } else if (job == 4) {
            mumps_S(3);
        }
    }while(true);

}		/* -----  end of function abcd::waitForSolve()  ----- */
