#include <abcd.h>
#include <iostream>
#include <fstream>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::solveS
 *  Description:  
 * =====================================================================================
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

        buildS(S_rows, S_cols, S_vals);

        if(inter_comm.rank() == 0){
            LDEBUG << "> T.build S : (no stop)  " << setprecision(2) << MPI_Wtime() - t;
        }

        // debug right now:
        bool write_sub_s = false;
        if(write_sub_s) {
            ofstream f;
            ostringstream ff;
            LINFO << "> Writing to file '" << "/tmp/m" << IRANK << ".mtx' ";
            ff << "/tmp/m";
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

    
        /*-----------------------------------------------------------------------------
         *  MUMPS part
         *-----------------------------------------------------------------------------*/
        mumps_S.sym = 1;
        mumps_S.par = 1;

        int job = 2;
        mpi::broadcast(intra_comm, job, 0);
        mumps_S.comm_fortran = MPI_Comm_c2f((MPI_Comm) comm);

        mumps_S(-1);
        mumps_S.initialized = true;

        mumps_S.icntl[0] = -1;
        mumps_S.icntl[1] = -1;
        mumps_S.icntl[2] = -1;

        if(inter_comm.rank() == 0){ 
            if(write_s.length() != 0)
                strcpy(mumps_S.write_problem, write_s.c_str());
            if(verbose == 2){
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

        t = MPI_Wtime();

        mumps_S(1);

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

        t = MPI_Wtime();

        mumps_S(2);

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
        /*-----------------------------------------------------------------------------
         *  END MUMPS part
         *-----------------------------------------------------------------------------*/
    } else {
        // tell the workers that we are going to solve with S, no building necessary!
        int job = 4;
        mpi::broadcast(intra_comm, job, 0);
    }
    

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
