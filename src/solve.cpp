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

    double t, tto = MPI_Wtime();
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
        mpi::all_reduce(inter_comm, f_ptr, size_c, f_o, or_bin);
        f = ff;
    }
    if(inter_comm.rank() == 0) cout << "Time to centralize f : " << MPI_Wtime() - t << endl;

    if(inter_comm.rank() == 0){
        cout << "*----------------------------------*" << endl;
        cout << "| [->] Solving Sz = f              |" << endl;
        cout << "*                                  *" << endl;
    }

    t = MPI_Wtime();
    if(icntl[15] != 0){
        if(inter_comm.rank() == 0)
            cout << "* ITERATIVELY                      *" << endl;

        MV_ColMat_double Ytf(m_l, 1, 0);

        double *f_ptr = f.ptr();

        VECTOR_double f0 = f.data();
        f0 = pcgS(f0);
        f.setData(f0);

        //mpi::broadcast(inter_comm, f_ptr, size_c, 0);

    } else {
        f = solveS(f);
    }
    if(IRANK == 0) 
        cout << "| Time to solve Sz = f : " << MPI_Wtime() - t << endl;

    t = MPI_Wtime();
    if(inter_comm.rank() == 0){
        cout << "*----------------------------------*" << endl;
        cout << "|                               T  |" << endl;
        cout << "| [->] Computing zz = (I - P) Y  z |" << endl;
        cout << "*----------------------------------*" << endl;
    }


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
    if(IRANK == 0) 
        cout << "| Other stuffs : " << MPI_Wtime() - t << endl;

    use_xk = false;

    // the final solution (distributed)
    f = w + f;
    if(IRANK == 0) cout << "Total time to build and solve " << MPI_Wtime() - tto << endl;

    double rho = compute_rho(f, b, 0);
    if(IRANK == 0) cout << "rho = " << rho << endl;

    //cout << "last element of f : " << f(n-1, 0) << endl;
    //if(IRANK == 0) cout << f(MV_VecIndex(0, 10), 0) << endl;

     //centralize the solution to the master
    {

        MV_ColMat_double xfmf(n, 1, 0);
        MV_ColMat_double lf(n, 1, 0);

        for(std::map<int,int>::iterator it = glob_to_local.begin(); it != glob_to_local.end(); it++){

            if(it->first < n_o){
                xfmf(it->second , 0) = 1 - f(it->second, 0);
                lf(it->second , 0) = f(it->second, 0);
            }
        }
        cout << "f: " <<  infNorm(lf) << endl;

        double nf = infNorm(xfmf);
        double nff = infNorm(lf);
        
        double nfa, nfb;
        mpi::all_reduce(inter_comm, &nf, 1, &nfa, mpi::maximum<double>());

        if(IRANK == 0) cout << "fwd : " <<  nfa << endl;
    }

    //if(IRANK == 0) cout << std::setprecision(16) << f << endl;
    //double rho;
    if(inter_comm.rank()==0){
        zrhs = MV_ColMat_double(m, 1, 0);
        int st = 0;
        for(int p = 0; p < partitions.size(); p++){
            zrhs(MV_VecIndex(st, st + partitions[p].dim(0) - 1),
                    MV_VecIndex(0, 0)) = spsmv(partitions[p], local_column_index[p], f);
            st += partitions[p].dim(0);
        }
        zrhs = zrhs - b;
        cout << "||\\bar{A}z - b||_2 = " << sqrt(zrhs.squaredSum()) << endl;
    }



}		/* -----  end of function abcd::solveABCD  ----- */


