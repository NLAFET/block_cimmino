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
 * \file sutils/solveS.cpp
 * \brief Build S and compute z=S^{-1}f using MUMPS
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include <abcd.h>
#include <iostream>
#include <fstream>

/*!
 *  \brief Build S and compute z=S^{-1}f using MUMPS
 *
 *  Build S and compute z=S^{-1}f using MUMPS
 *
 *  \param f: the right hand side
 *
 */
MV_ColMat_double
abcd::solveS ( MV_ColMat_double &f )
{
    double t;


    if(inter_comm.rank() == 0){
        LINFO << "*----------------------------------*";
        LINFO << "> Building S = Y (I - P) Y^T       |";
        LINFO << "*----------------------------------*";
    }

    // if not created yet, do it!
    if( !mumps_S.initialized ){
        t = MPI_Wtime();

        // build complete S matrix
        buildS(S_rows, S_cols, S_vals);

        if(inter_comm.rank() == 0){
            LDEBUG << "> T.build S : (no stop)  " << setprecision(2) << MPI_Wtime() - t;
        }

#ifdef WIP
        // debug right now:
        if(write_s.length() != 0) {
            ofstream f;
            ostringstream ff;
            LINFO << "> Writing to file '" << write_s << "_" << IRANK << ".mtx' ";
            ff << write_s << "_";
            ff << IRANK << ".mtx";
            f.open(ff.str().c_str());
            f << "%%MatrixMarket matrix coordinate real general\n";
            f << size_c << " " << size_c << " " << S_rows.size() << "\n";
            for(size_t i = 0; i < S_vals.size(); i++){
                f << S_rows[i] << " " << S_cols[i] << " " << S_vals[i] << "\n";
            }
            f.close();
            exit(0);
        }
#endif //WIP

        /*-----------------------------------------------------------------------------
         *  MUMPS initialization
         *-----------------------------------------------------------------------------*/
        mumps_S.sym = 1;
        mumps_S.par = 1;

        int job = 2;
        mpi::broadcast(intra_comm, job, 0);
        mumps_S.comm_fortran = MPI_Comm_c2f((MPI_Comm) comm);

        mumps_S(-1);
        mumps_S.initialized = true;

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

        if(inter_comm.rank() == 0){
            if(write_s.length() != 0)
                strcpy(mumps_S.write_problem, write_s.c_str());
            if(icntl[Controls::mumps_verbose]){
                mumps_S.icntl[0] = 6;
                mumps_S.icntl[1] = 6;
                mumps_S.icntl[2] = 6;
            }
        }

        mumps_S.n = size_c;

        mumps_S.setIcntl(8, 77);
        mumps_S.setIcntl(7, 5);
        mumps_S.setIcntl(14, 70);
        mumps_S.setIcntl(40, 1);

        // Initialize matrix for MUMPS
        if(inter_comm.size() == 1){
            mumps_S.nz= S_nnz();
            mumps_S.irn = &S_rows[0];
            mumps_S.jcn = &S_cols[0];
            mumps_S.a= &S_vals[0];

        } else {
            mumps_S.icntl[18 - 1]= 3;
            mumps_S.nz_loc = S_nnz();

            mumps_S.irn_loc = &S_rows[0];
            mumps_S.jcn_loc = &S_cols[0];
            mumps_S.a_loc   = &S_vals[0];
        }

        if(inter_comm.rank() == 0){
            LINFO << "*                                  *";
            LINFO << "> T. build and init mumps S: " << setprecision(2) << MPI_Wtime() - t;
            LINFO << "*                                  *";
        }

        /*-----------------------------------------------------------------------------
         *  MUMPS analysis
         *-----------------------------------------------------------------------------*/
        t = MPI_Wtime();

        // launch MUMPS analysis
        mumps_S(1);

        // check MUMPS analysis
        if(mumps_S.info[0] < 0) {
            LERROR << "MUMPS exited with " << mumps_S.info[0];
            int job = -90 + mumps_S.info[0];
            mpi::broadcast(intra_comm, job, 0);
            throw std::runtime_error("MUMPS exited with an error");
        }

        if(inter_comm.rank() == 0){
            LINFO << "*                                  *";
            LINFO << "> T.Analyse S:   " << setprecision(2) << MPI_Wtime() - t;
            LINFO << "*                                  *";
        }

        /*-----------------------------------------------------------------------------
         *  MUMPS factorization
         *-----------------------------------------------------------------------------*/
        t = MPI_Wtime();

        // launch MUMPS facto
        mumps_S(2);

        // check MUMPS facto
        if(mumps_S.info[0] < 0) {
            LERROR << "MUMPS exited with " << mumps_S.info[0];
            int job = -90 + mumps_S.info[0];
            mpi::broadcast(intra_comm, job, 0);
            throw std::runtime_error("MUMPS exited with an error");
        }

        if(inter_comm.rank() == 0){
            LINFO << "*                                  *";
            LINFO << "> T.Factorize S: " << setprecision(2) << MPI_Wtime() - t;
            LINFO << "*                                  *";
        }
    } else {
        // tell the workers that we are going to solve with S, no building necessary!
        //todo: this should only happen if MUMPS didn't factorize the matrix!
        int job = 4;
        mpi::broadcast(intra_comm, job, 0);
    }


    /*-----------------------------------------------------------------------------
     *  MUMPS Solve
     *-----------------------------------------------------------------------------*/
    t = MPI_Wtime();

    if(inter_comm.rank() == 0){
        mumps_S.rhs = f.ptr();
        mumps_S.nrhs = 1;
        mumps_S.lrhs = size_c;
    }

    if(inter_comm.rank() == 0){
        LINFO << "*                                  *";
        LINFO << "> Computing z = S^-1 f             |";
    }

    mumps_S(3);

    double *f_ptr = f.ptr();
    // TODO : better send parts not the whole z
    //
    mpi::broadcast(inter_comm, f_ptr, size_c, 0);
    if(inter_comm.rank() == 0){
        LINFO << "> Took: " << setprecision(2) << MPI_Wtime() - t;
        LINFO << "*----------------------------------*";
    }

    return f;
}       /* -----  end of function abcd::solveS  ----- */
