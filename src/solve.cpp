/*
 * =====================================================================================
 *
 *       Filename:  solve.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/06/2013 11:36:33 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */


#include <abcd.h>
#include <iostream>
#include <fstream>


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::solveABCD
 *  Description:  
 * =====================================================================================
 */
    void
abcd::solveABCD ( MV_ColMat_double &b )
{

    double t;
    MV_ColMat_double w;
    if(inter_comm.rank() == 0){
        cout << "*----------------------------------*" << endl;
        cout << "| [->] Computing w = A^+b          |" << endl;
        cout << "*----------------------------------*" << endl;
    }

    t = MPI_Wtime();
    if(dcntl[10] == 0){
        w = sumProject(1e0, b, 0e0, Xk);
    } else {
        bcg(b);
        w = Xk; 
    }
    if(inter_comm.rank() == 0) cout << "Time to compute w = A^+ b : " << MPI_Wtime() - t << endl;

    if(inter_comm.rank() == 0){
        cout << "*----------------------------------*" << endl;
        cout << "| [->] Building S = Y (I - P) Y^T  |" << endl;
        cout << "*----------------------------------*" << endl;
    }

    // if not created yet, do it!
    t = MPI_Wtime();
    if( S.dim(0) == 0 ){
        S = abcd::buildS();
    }

    inter_comm.barrier(); //useless! used for timing
    if(inter_comm.rank() == 0) cout << "Time to build S : " << MPI_Wtime() - t << endl;

    //if(inter_comm.rank() == 0) cout << S(0,0) << endl;
    //inter_comm.barrier();
    //
    
    /*-----------------------------------------------------------------------------
     *  MUMPS part
     *-----------------------------------------------------------------------------*/
    t = MPI_Wtime();
    DMUMPS_STRUC_C mu;
    mu.sym = 0;
    mu.par = 1;
    mu.job = -1;
    mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) inter_comm);

    dmumps_c(&mu);

    mu.icntl[0] = -1;
    mu.icntl[1] = -1;
    mu.icntl[2] = -1;

    //if(inter_comm.rank() == 0){ 
        //strcpy(mu.write_problem, "/tmp/chiante.mtx");
        //mu.icntl[0] = 6;
        //mu.icntl[1] = 6;
        //mu.icntl[2] = 6;
        //mu.icntl[3] = 2;
    //}

    mu.n = S.dim(0);

    // parallel analysis if the S is large enough
    //if(mu.n >= 200) {
        //mu.icntl[28 - 1] =  2;
    //}
    mu.icntl[8  - 1] =  7;
    mu.icntl[7  - 1] =  5;
    mu.icntl[14 - 1] =  70;

    if(inter_comm.size() == 1){ 
        mu.nz= S.NumNonzeros();
        //mu.irn= S.rowind_ptr();
        //mu.jcn= S.colind_ptr();
        mu.irn = new int[mu.nz];
        mu.jcn = new int[mu.nz];
        mu.a= S.val_ptr();
        for(int i = 0; i < mu.nz; i++){
            //mu.irn[i]++;
            //mu.jcn[i]++;
            mu.irn[i] = S.row_ind(i) + 1;
            mu.jcn[i] = S.col_ind(i) + 1;
        }
    } else {
        mu.icntl[18 - 1]= 3;
        mu.nz_loc = S.NumNonzeros();

        mu.irn_loc = S.rowind_ptr();
        mu.jcn_loc = S.colind_ptr();
        mu.a_loc = S.val_ptr();
        for(int i = 0; i < mu.nz_loc; i++){
            mu.irn_loc[i]++;
            mu.jcn_loc[i]++;
        }
    }
    t = MPI_Wtime();

    mu.job = 1;
    dmumps_c(&mu);

    if(inter_comm.rank() == 0) cout << "Time to Analyse S : " << MPI_Wtime() - t << endl;

    t = MPI_Wtime();

    mu.job = 2;
    dmumps_c(&mu);

    if(inter_comm.rank() == 0) cout << "Time to Factorize S : " << MPI_Wtime() - t << endl;

    if(mu.info[0] < 0) {
        cout << mu.info[0] << endl;
        exit(0); 
    }
    /*-----------------------------------------------------------------------------
     *  END MUMPS part
     *-----------------------------------------------------------------------------*/

    if(inter_comm.rank() == 0){
        cout << "*----------------------------------*" << endl;
        cout << "| [->] Setting f = - Y w           |" << endl;
        cout << "*----------------------------------*" << endl;
    }

    t = MPI_Wtime();
    MV_ColMat_double f(size_c, 1, 0);
    for(std::map<int,int>::iterator it = glob_to_local.begin(); it != glob_to_local.end(); it++){
        if(it->first >= n_o){
            f(it->first - n_o, 0) = -1 * w(it->second, 0);
        }
    }
    {
        double *f_ptr = f.ptr();
        MV_ColMat_double ff(size_c, 1, 0);
        double *f_o = ff.ptr();
        mpi::reduce(inter_comm, f_ptr, size_c, f_o, or_bin, 0);
        f = ff;
    }
    if(inter_comm.rank() == 0) cout << "Time to centralize f : " << MPI_Wtime() - t << endl;

    if(inter_comm.rank() == 0){
        mu.rhs = f.ptr();
        mu.nrhs = 1;
    }

    if(inter_comm.rank() == 0){
        cout << "*----------------------------------*" << endl;
        cout << "| [->] Solving Sz = f              |" << endl;
        cout << "*----------------------------------*" << endl;
    }

    mu.job = 3;
    dmumps_c(&mu);

    //if(inter_comm.size() == 1){ 
        //for(int i = 0; i < mu.nz; i++){
            //mu.irn[i]--;
            //mu.jcn[i]--;
        //}
    //} else {
        //for(int i = 0; i < mu.nz_loc; i++){
            //mu.irn_loc[i]--;
            //mu.jcn_loc[i]--;
        //}
    //}
    //

    if(inter_comm.rank() == 0){
        cout << "*----------------------------------*" << endl;
        cout << "|                               T  |" << endl;
        cout << "| [->] Computing zz = (I - P) Y  z |" << endl;
        cout << "*----------------------------------*" << endl;
    }


    // broadcast f to other cpus, where f is the new z
    double *f_ptr = f.ptr();
    // TODO : better send parts not the whole z
    //
    mpi::broadcast(inter_comm, f_ptr, size_c, 0);


    Xk = MV_ColMat_double(n, 1, 0);
    MV_ColMat_double zrhs(m, 1, 0); 

    for( int i = 0; i < size_c; i++){

        std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);

        if(iti!=glob_to_local.end()){
            Xk(glob_to_local[n_o + i], 0) = f(i, 0);
        } else {
            continue;
        }
    }

    f = MV_ColMat_double(n, 1, 0);

    if(dcntl[10] == 0){
        f = Xk - sumProject(0e0, b, 1e0, Xk);

    } else {
        use_xk = false;
        //itmax = 2;

        if(!use_xk){
            int st = 0;
            for(int p = 0; p < partitions.size(); p++){
                //zrhs(MV_VecIndex(st, st + partitions[p].dim(0) - 1),
                        //MV_VecIndex(0, 0)) = spsmv(partitions[p], local_column_index[p], Xk);
                MV_ColMat_double sp =  spsmv(partitions[p], local_column_index[p], Xk);
                int pos = 0;
                for(int k = st; k < st + partitions[p].dim(0); k++){
                    zrhs(k, 0) = sp(pos, 0);
                    pos++;
                }
                st += partitions[p].dim(0);
            }

            for(int i = 0; i < f.dim(0); i++)
                f(i, 0) = Xk(i, 0);
        }

        bcg(zrhs);

        if(use_xk)
            f = Xk;
        else
            f = f - Xk;
    }

    use_xk = false;

    // the final solution (distributed)
    f = w + f;

    double rho = compute_rho(f, b, 0);
    if(IRANK == 0) cout << "rho = " << rho << endl;

    //if(inter_comm.rank()==0){
        //zrhs = MV_ColMat_double(m, 1, 0);
        //int st = 0;
        //for(int p = 0; p < partitions.size(); p++){
            //zrhs(MV_VecIndex(st, st + partitions[p].dim(0) - 1),
                    //MV_VecIndex(0, 0)) = spsmv(partitions[p], local_column_index[p], f);
            //st += partitions[p].dim(0);
        //}
        //zrhs = zrhs - b;
        //cout << "||Ax - b||_2 = " << sqrt(zrhs.squaredSum()) << endl;
    //}


    // centralize the solution to the master
    //{
        //double *f_ptr = f.ptr();
        //MV_ColMat_double ff(size_c, 1, 0);
        //double *f_o = ff.ptr();
        //mpi::reduce(inter_comm, f_ptr, size_c, f_o, or_bin, 0);
        //f = ff;
    //}




    //// set xk = [0 I] 
    //Xk = MV_ColMat_double(n, 1, 0);
    //// Retrieve all keys
    ////cout << inter_comm.rank() << " " << indices[0] << " " << indices.back() << " " << n_o + i << endl;
    //std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);

    //if(iti!=glob_to_local.end()){
        //Xk(glob_to_local[n_o + i], 0) = 1;
    //}

    //MV_ColMat_double b(m, 1, 0); 

}		/* -----  end of function abcd::solveABCD  ----- */


