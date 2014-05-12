#include <abcd.h>
#include <iostream>
#include <fstream>

void abcd::solveABCD ( MV_ColMat_double &b )
{
    {
        MV_ColMat_double u(m, nrhs, 0);
        u = b(MV_VecIndex(0, b.dim(0)-1), MV_VecIndex(0,nrhs-1));
        double lnrmBs = infNorm(u);
        // Sync B norm :
        mpi::all_reduce(inter_comm, &lnrmBs, 1,  &nrmB, mpi::maximum<double>());
        mpi::broadcast(inter_comm, nrmMtx, 0);
    }

    double t, tto = MPI_Wtime();
    MV_ColMat_double w;
    if(inter_comm.rank() == 0){
        LINFO << "*----------------------------------*";
        LINFO << "> Computing w = A^+b                ";
    }

    t = MPI_Wtime();
    if(dcntl[Controls::aug_filter] == 0){
        w = sumProject(1e0, b, 0e0, Xk);
    } else {
        bcg(b);
        w = Xk; 
    }
    if(inter_comm.rank() == 0){
        LINFO << "> Time to compute w = A^+ b: " << setprecision(2) << MPI_Wtime() - t;
    }
    

    if(inter_comm.rank() == 0){
        LINFO << "*----------------------------------*";
        LINFO << "> Setting f = - Y w                 ";
    }

    t = MPI_Wtime();
    MV_ColMat_double f(size_c, 1, 0);
    for(std::map<int,int>::iterator it = glob_to_local.begin();
        it != glob_to_local.end(); ++it){
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
    if(inter_comm.rank() == 0){
        LINFO << "> Time to centralize f : " << setprecision(2) << MPI_Wtime() - t;
    }


    if(inter_comm.rank() == 0){
        LINFO << "*----------------------------------*";
        LINFO << "> Solving Sz = f                    ";
        LINFO << "*                                  *";
    }

    t = MPI_Wtime();
    if(icntl[Controls::aug_iterative] != 0){
        if(inter_comm.rank() == 0)
            LINFO << "* ITERATIVELY                      *";

        VECTOR_double f0 = f.data();
        f0 = pcgS(f0);
        f.setData(f0);

    } else {

        f = solveS(f);
    }
    if(IRANK == 0) 
        LINFO << "| Time to solve Sz = f: " << setprecision(2) << MPI_Wtime() - t;

    t = MPI_Wtime();
    if(inter_comm.rank() == 0){
        LINFO << "*----------------------------------*";
        LINFO << "                                T   ";
        LINFO << "> [->] Computing zz = (I - P) Y  z  ";
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

    if(dcntl[Controls::aug_filter] == 0){
        f = Xk - sumProject(0e0, b, 1e0, Xk);

    } else {
        use_xk = false;

        if(!use_xk){
            int st = 0;
            for(int p = 0; p < nb_local_parts; p++){
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
        LINFO << "> Took: " << setprecision(2) << MPI_Wtime() - t;


    //! the final solution (distributed)
    //! w = \Abar^+ b
    //! f = \Wbar^+ b
    Xk = w + f;

    if(IRANK == 0){
        LINFO << "*----------------------------------*";
        LINFO << "> Total time to solve: " << setprecision(2) << MPI_Wtime() - tto;
        LINFO << "*----------------------------------*";
        LINFO << "";
    }

    compute_rho(Xk, b);

    if(IRANK == 0) {
        t = MPI_Wtime();
        LINFO << "> Centralizing solution";
        MV_ColMat_double temp_sol(n_o, 1, 0);
        std::map<int, std::vector<double> > xo;
        std::map<int, std::vector<int> > io;
        for(int k = 1; k < inter_comm.size(); k++){
            inter_comm.recv(k, 71, xo[k]);
            inter_comm.recv(k, 72, io[k]);
        }
        for(int k = 1; k < inter_comm.size(); k++){
            for(size_t i = 0; i < io[k].size() && io[k][i] < n_o; i++){
                int ci = io[k][i];
                temp_sol(ci, 0) = xo[k][i] * dcol_(io[k][i]);
            }
        }

        for(size_t i = 0; i < glob_to_local_ind.size() && glob_to_local_ind[i] < n_o ; i++){
            temp_sol(glob_to_local_ind[i], 0) = Xk(i, 0) * dcol_(glob_to_local_ind[i]);
        }

        sol = temp_sol(MV_VecIndex(0, n_o-1), 0);

        if(Xf.dim(0) != 0) {
            MV_ColMat_double xf = MV_ColMat_double(n, 1, 0);
            xf = Xf - sol;
            double nrmxf =  infNorm(xf);
            dinfo[Controls::forward_error] = nrmxf/nrmXf;
        }
        LINFO << "> Took: " << setprecision(2) << MPI_Wtime() - t;

    } else {
        std::vector<double> x;
        x.reserve(n);
        for(int i = 0; i < n; i++){
            x.push_back(Xk(i, 0));
        }
        inter_comm.send(0, 71, x);
        inter_comm.send(0, 72, glob_to_local_ind);
    }

}		/* -----  end of function abcd::solveABCD  ----- */


