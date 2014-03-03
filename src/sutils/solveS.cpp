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

    // debug right now:
    bool write_sub_s = false;

    if(inter_comm.rank() == 0){
        clog << "*      ----------------------      *" << endl;
        clog << "| [--] Building S = Y (I - P) Y^T  |" << endl;
        clog << "*      ----------------------      *" << endl;
    }

    // if not created yet, do it!
    t = MPI_Wtime();
    if( S_vals.size() == 0 ){
        buildS(S_rows, S_cols, S_vals);
    }

    if(inter_comm.rank() == 0){
        clog << "*                                  *" << endl;
        clog << "  [--] T.build S : (no stop)  " << MPI_Wtime() - t << endl;
        clog << "*                                  *" << endl;
    }

    if(write_sub_s) {
        ofstream f;
        ostringstream ff;
        ff << "/tmp/m";
        ff << IRANK << ".mtx";
        f.open(ff.str().c_str());
        f << "%%MatrixMarket matrix coordinate real general\n";
        f << S.dim(0) << " " << S.dim(1) << " " << S.NumNonzeros() << "\n";
        for(int i = 0; i < S.NumNonzeros(); i++){
            f << S.row_ind(i) + 1 << " " << S.col_ind(i) + 1 << " " << S.val(i) << "\n";
        }
        f.close();
    }

    
    /*-----------------------------------------------------------------------------
     *  MUMPS part
     *-----------------------------------------------------------------------------*/
    MUMPS mu;
    mpi::communicator world;
    mu.sym = 2;
    mu.par = 1;
    mu.job = -1;

    int job = 2;
    mpi::broadcast(intra_comm, job, 0);
    mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) world);

    dmumps_c(&mu);

    mu.icntl[0] = -1;
    mu.icntl[1] = -1;
    mu.icntl[2] = -1;

    if(inter_comm.rank() == 0){ 
        if(write_s.length() != 0)
            strcpy(mu.write_problem, write_s.c_str());
        mu.icntl[0] = -1;
        mu.icntl[1] = -1;
        mu.icntl[2] = -1;
    }

    mu.n = size_c;

    mu.setIcntl(8, 77);
    mu.setIcntl(7, 5);
    mu.setIcntl(14, 70);
    mu.setIcntl(40, 1);

    if(inter_comm.size() == 1){ 
        mu.nz= S_nnz();
        mu.irn = &S_rows[0];
        mu.jcn = &S_cols[0];
        mu.a= &S_vals[0];

    } else {
        mu.icntl[18 - 1]= 3;
        mu.nz_loc = S_nnz();

        mu.irn_loc = &S_rows[0];
        mu.jcn_loc = &S_cols[0];
        mu.a_loc   = &S_vals[0];
    }

    if(inter_comm.rank() == 0){
        clog << "*                                  *" << endl;
        clog << "  [--] T build and init mumps S: " << MPI_Wtime() - t << endl;
        clog << "*                                  *" << endl;
    }

    t = MPI_Wtime();

    mu.job = 1;
    dmumps_c(&mu);

    if(inter_comm.rank() == 0){
        clog << "*                                  *" << endl;
        clog << "  [--] T.Analyse S :   " << MPI_Wtime() - t << endl;
        clog << "*                                  *" << endl;
    }

    t = MPI_Wtime();

    mu.job = 2;
    dmumps_c(&mu);

    if(inter_comm.rank() == 0){
        clog << "*                                  *" << endl;
        clog << "  [--] T.Factorize S : " << MPI_Wtime() - t << endl;
        clog << "*                                  *" << endl;
    }

    if(mu.info[0] < 0) {
        clog << mu.info[0] << endl;
        exit(0); 
    }
    /*-----------------------------------------------------------------------------
     *  END MUMPS part
     *-----------------------------------------------------------------------------*/


    if(inter_comm.rank() == 0){
        mu.rhs = f.ptr();
        mu.nrhs = 1;
    }

    if(inter_comm.rank() == 0){
        clog << "*      -------------------         *" << endl;
        clog << "| [--] Computing z = S^-1 f        |" << endl;
        clog << "*----------------------------------*" << endl;
    }

    mu.job = 3;
    dmumps_c(&mu);

    double *f_ptr = f.ptr();
    // TODO : better send parts not the whole z
    //
    mpi::broadcast(inter_comm, f_ptr, size_c, 0);

    return f;
}       /* -----  end of function abcd::solveS  ----- */
