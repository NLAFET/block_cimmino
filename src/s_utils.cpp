/*
 * =====================================================================================
 *
 *       Filename:  s_utils.cpp
 *
 *    Description:  Contain the different methods related to building S and
 *    computing columns of S
 *
 *        Version:  1.0
 *        Created:  03/11/2013 02:11:05 PM
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
 *         Name:  abcd::buildS
 *  Description:  Builds a sparse S
 * =====================================================================================
 */
    Coord_Mat_double
abcd::buildS (  )
{
    std::vector<int> cols;
    return buildS( cols );
}
    Coord_Mat_double
abcd::buildS ( std::vector<int> cols )
{
    Coord_Mat_double shur;
    //MV_ColMat_double S(size_c, size_c, 0);
    std::vector<int> sr, sc;
    std::vector<double> sv;

    std::vector<int> vc, vr;
    std::vector<double> vv;

    std::vector<int> my_cols;

    if(cols.size() == 0) {
        for( int i = 0; i < size_c; i++){
            std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);
            if(iti!=glob_to_local.end()) my_cols.push_back(i);
        }
    } else {
        for( size_t i = 0; i < cols.size(); i++){
            std::map<int,int>::iterator iti = glob_to_local.find(n_o + cols[i]);
            if(iti!=glob_to_local.end()) my_cols.push_back(cols[i]);
        }
    }

    int maxcols = mpi::all_reduce(inter_comm, (int) my_cols.size(), mpi::maximum<int>());
    int mincols = mpi::all_reduce(inter_comm, (int) my_cols.size(), mpi::minimum<int>());
    int total = mpi::all_reduce(inter_comm, (int) my_cols.size(), std::plus<int>());
    IFMASTER clog << "Max number of cols is " << maxcols << " ,  min is " << mincols <<
       "and average is " << total/parallel_cg <<    endl;

#ifdef MUMPS_ES
    mumps.keep[235 - 1] = icntl[8];
    mumps.keep[261 - 1] = icntl[8];
    mumps.keep[495 - 1] = icntl[8];
    mumps.keep[497 - 1] = icntl[8];
#endif
    if(dcntl[10] == 0 || icntl[15] == 2){

        std::vector<int>::iterator pos = my_cols.begin();
        std::vector<int>::iterator end_pos;

        vc.reserve(my_cols.size() * my_cols.size());
        vr.reserve(my_cols.size() * my_cols.size());
        vv.reserve(my_cols.size() * my_cols.size());

        int share = icntl[14];
        while(pos != my_cols.end()){
            if(pos + share < my_cols.end()) end_pos = pos + share;
            else end_pos = my_cols.end();
            
            double perc = end_pos - my_cols.begin();
            perc /= my_cols.size();
            perc *= 100;

            //clog << IRANK << " ["<< floor(perc) << "%] : " 
                //<< end_pos - my_cols.begin() << " / " << my_cols.size()
                //<< endl;

            std::vector<int> cur_cols;

            std::copy(pos, end_pos, std::back_inserter(cur_cols));


            //int mumps_share = share > 32 ? share : 16;
            int mumps_share = share;
            mumps.icntl[27 - 1] = mumps_share;

            if(intra_comm.size() == 1){
                MV_ColMat_double sp = spSimpleProject(cur_cols);

                double *sptr = sp.ptr();
                int slda = sp.lda();
                int srows = sp.dim(0);

                for( size_t j = 0; j < cur_cols.size(); j++){
                    int c = cur_cols[j];
                    for( int i = c; i < srows; i++){
                        if(sptr[i + j * slda] != 0){
                        //if(sp(i,j)!=0){
                        //if(abs(sp(i,j))>1e-16){
                            vc.push_back(c);
                            vr.push_back(i);
                            vv.push_back(sptr[i + j * slda]);
                        }
                    }
                }
            } else {
                vector<int> target_idx;
                vector<int> target;
                target.reserve(size_c);
                target_idx.reserve(size_c);
                spSimpleProject(cur_cols, vr, vc, vv, target, target_idx);
            }
            pos = end_pos;
        }


    } else {

        //for( int j = 0; j < my_cols.size(); j++){
            //int i = my_cols[j];
        for( int i = 0; i < size_c; i++){
        //
            //int my_bro = -1;
            if(inter_comm.rank() == 0) clog << " Column " << i << " out of " << size_c << endl;

            //for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
                    //it != col_interconnections.end(); it++) {
                //int k = it->second.size() - 1;
                //while( k > -1 && my_bro == -1){
                    //if(it->second[k] == glob_to_local[n_o + i]){
                        //my_bro = it->first;
                    //}
                    //k--;
                //}
            //}


            // set xk = [0 I] 
            block_size = 1;
            Xk = MV_ColMat_double(n, 1, 0);
            MV_ColMat_double b(m, 1, 0); 

            std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);
            if(iti!=glob_to_local.end()){
                Xk(glob_to_local[n_o + i], 0) = 1;

                //int st = 0;
                //for(int p = 0; p < nb_local_parts; p++){
                    //MV_ColMat_double sp = spsmv(partitions[p], local_column_index[p], Xk);
                    //int pos = 0;
                    //for(int k = st; k < st + partitions[p].dim(0); k++){
                        //b(k, 0) = sp(pos, 0);
                        //pos++;
                    //}
                    //st += partitions[p].dim(0);
                //}
            }
            //MV_ColMat_double sp = Xk - coupleSumProject(0e0, b, 1e0, Xk, my_bro);

            //itmax = 3;
            use_xk = true;
            //verbose = true;
            //threshold = 1e-8;
            bcg(b);
            use_xk = false;

            for(std::map<int,int>::iterator it = glob_to_local.begin(); it != glob_to_local.end(); it++){

                if(it->first >= n_o && comm_map[it->second] == 1){
                    //if(inter_comm.rank() > my_bro) continue;

                    //clog << "Comm map : " << it->first << " " << comm_map[it->first] << endl;
                    vc.push_back(i);
                    vr.push_back(it->first - n_o);
                    vv.push_back(Xk(it->second,0));
                }
            }
        }
    }
    if(vv.size() == 0) {
        vc.push_back(0);
        vr.push_back(0);
        vv.push_back(0.0);
    }

    shur = Coord_Mat_double(size_c, size_c, vv.size(), &vv[0], &vr[0], &vc[0]);

    use_xk = false;

    return shur;
}       /* -----  end of function abcd::buildS  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::prodSv
 *  Description:  
 * =====================================================================================
 */
    MV_ColMat_double
abcd::prodSv ( MV_ColMat_double &V )
{
    MV_ColMat_double W(n, V.dim(1), 0);
    MV_ColMat_double b; // just for decoration!
    MV_ColMat_double R(V.dim(0), V.dim(1), 0);

    for(std::vector<int>::iterator iti = st_c_part_it; iti != glob_to_local_ind.end(); iti++){
            for(int j = 0; j < V.dim(1); j++)
                W(iti - glob_to_local_ind.begin() , j) = V(*iti - n_o, j);
    }
    //for(std::map<int,int>::iterator it = glob_to_local.begin(); it != glob_to_local.end(); it++){
        //if(it->first >= n_o){
                //for(int j = 0; j < V.dim(1); j++)
                    //W(it->second , j) = V(it->first - n_o, j);
        //}
    //}

    W = W - sumProject(0e0, b, 1e0, W);
    //for(std::map<int,int>::iterator it = glob_to_local.begin(); it != glob_to_local.end(); it++){
        //if(it->first >= n_o){
                    //for(int j = 0; j < V.dim(1); j++) R(it->first - n_o, j) = W(it->second , j); 
        //}
    //}

    for(std::vector<int>::iterator iti = st_c_part_it; iti != glob_to_local_ind.end(); iti++){
            for(int j = 0; j < V.dim(1); j++)
                R(*iti - n_o, j) = W(iti - glob_to_local_ind.begin() , j);

    }

    return R;
}       /* -----  end of function abcd::prodSv  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::buildM
 *  Description:  
 * =====================================================================================
 */
MUMPS abcd::buildM (  )
{
    if(inter_comm.rank() == 0){
        clog << "* Building the preconditioner      *" << endl;
        clog << " Size is : " << selected_S_columns.size() << endl;

        //ofstream f;
        //f.open("/tmp/selected");
        //for(int i = 0; i < selected_S_columns.size(); i++){
            //f << selected_S_columns[i] << endl;
        //}
        //f.close();
    }


    double t = MPI_Wtime();
    //t = MPI_Wtime();
    //S = buildS(skipped_S_columns);
    //clog << " Time to build S : " << MPI_Wtime() - t << endl;

    Coord_Mat_double M = buildS(selected_S_columns);
    //clog << " Time to build M : " << MPI_Wtime() - t << endl;

    //if(IRANK == 0) clog << "* -----------  DONE  ------------- * " << endl;

    std::vector<int> dropped;

    {
        std::vector<int> mr, mc;
        std::vector<double> mv;

        std::vector<int>::iterator it;

        for(size_t i = 0; i < skipped_S_columns.size(); i++){
            std::map<int,int>::iterator iti = glob_to_local.find(n_o + skipped_S_columns[i]);
            if(iti == glob_to_local.end()) continue;

            int ro = skipped_S_columns[i];

            mr.push_back(ro);
            mc.push_back(ro);
            //mv.push_back(S(ro,ro));

            mv.push_back(1);
        }

        for(int i = 0; i < M.NumNonzeros(); i++){
            if(std::find(skipped_S_columns.begin(), skipped_S_columns.end(), M.row_ind(i))
                    !=skipped_S_columns.end()) continue;
            if(M.val(i) == 0) continue;
            mr.push_back(M.row_ind(i));
            mc.push_back(M.col_ind(i));
            mv.push_back(M.val(i));
        }
        
        M = Coord_Mat_double(size_c, size_c, mv.size(), &mv[0], &mr[0], &mc[0]);
    }

    /*-----------------------------------------------------------------------------
     *  MUMPS part
     *-----------------------------------------------------------------------------*/
    
    MUMPS mu;
    mu.sym = 1;
    mu.par = 1;
    mu.job = -1;
    mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) inter_comm);

    dmumps_c(&mu);

    mu.icntl[0] = -1;
    mu.icntl[1] = -1;
    mu.icntl[2] = -1;

    if(inter_comm.rank() == 0){ 
        //strcpy(mu.write_problem, "/tmp/mmm.mtx");
        //mu.icntl[0] = 6;
        //mu.icntl[1] = 6;
        //mu.icntl[2] = 6;
        //mu.icntl[3] = 2;
    }

    mu.n = M.dim(0);

    // parallel analysis if the S is large enough
    //if(mu.n >= 200) {
        //mu.icntl[28 - 1] =  2;
    //}
    mu.icntl[8  - 1] =  7;
    mu.icntl[7  - 1] =  5;
    mu.icntl[14 - 1] =  90;

    if(inter_comm.size() == 1){ 
        mu.nz= M.NumNonzeros();
        mu.irn= M.rowind_ptr();
        mu.jcn= M.colind_ptr();
        mu.a= M.val_ptr();
        for(int i = 0; i < mu.nz; i++){
            mu.irn[i]++;
            mu.jcn[i]++;
        }
    } else {
        mu.icntl[18 - 1]= 3;
        mu.nz_loc = M.NumNonzeros();

        mu.irn_loc = M.rowind_ptr();
        mu.jcn_loc = M.colind_ptr();
        mu.a_loc = M.val_ptr();
        for(int i = 0; i < mu.nz_loc; i++){
            mu.irn_loc[i]++;
            mu.jcn_loc[i]++;
        }
    }
    t = MPI_Wtime();

    IBARRIER;
    if(inter_comm.rank() == 0)
        clog << "* Starting Analysis of the precond *" << endl;

    mu.job = 1;
    dmumps_c(&mu);

    if(inter_comm.rank() == 0){
        clog << "*                                  *" << endl;
        clog << "  [--] T.Analyse M :   " << MPI_Wtime() - t << endl;
        clog << "*                                  *" << endl;
    }

    t = MPI_Wtime();

    mu.job = 2;
    dmumps_c(&mu);

    if(inter_comm.rank() == 0){
        clog << "*                                  *" << endl;
        clog << "  [--] T.Factorize M : " << MPI_Wtime() - t << endl;
        clog << "  N, NZ              : " << mu.n << ", " << mu.nz << endl;
        clog << "*                                  *" << endl;
    }

    mu.icntl[0] = -1;
    mu.icntl[1] = -1;
    mu.icntl[2] = -1;

    if(mu.info[0] < 0) {
        clog << IRANK << " 1 : " << mu.info[0] <<  " 2 : " << mu.info[1] << endl;
        exit(0); 
    }
    /*-----------------------------------------------------------------------------
     *  END MUMPS part
     *-----------------------------------------------------------------------------*/
    return mu;
}       /* -----  end of function abcd::buildM  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::solveM
 *  Description:  
 * =====================================================================================
 */
    VECTOR_double
abcd::solveM (MUMPS &mu, VECTOR_double &z )
{
    if(IRANK == 0) {
        mu.rhs = new double[z.size()];
        for(int i = 0; i < z.size(); i++) mu.rhs[i] = z(i);
        mu.nrhs = 1;
    }

    mu.job = 3;
    dmumps_c(&mu);
    if(IRANK == 0){
        VECTOR_double sol(mu.rhs, z.size());
        double *sol_ptr = sol.ptr();
        mpi::broadcast(inter_comm, sol_ptr, size_c, 0);
        delete[] mu.rhs;
        return sol;
    } else {
        VECTOR_double sol(z.size(), 0);
        double *sol_ptr = sol.ptr();
        mpi::broadcast(inter_comm, sol_ptr, size_c, 0);
        return sol;
    }

}       /* -----  end of function abcd::solveM  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::pcgS
 *  Description:  
 * =====================================================================================
 */
    VECTOR_double
abcd::pcgS ( VECTOR_double &b )
{
    double resid, tol = 1e-5;

    //int max_iter = 2;
    int max_iter = size_c;

    MUMPS mu;
    if(dcntl[15] != 0) mu = buildM();
    double t = MPI_Wtime();

    VECTOR_double p, z, q;
    VECTOR_double alpha(1), beta(1), rho(1), rho_1(1);

    MV_ColMat_double zk(size_c, 1, 0);
    VECTOR_double x = zk.data();

    double normb = norm(b);
    MV_ColMat_double Szk = prodSv(zk);
    VECTOR_double r = b - Szk(0);

    if (normb == 0.0) 
        normb = 1;

    resid = norm(r) / normb;

    if (resid <= tol) {
        tol = resid;
        max_iter = 0;
        return 0;
    }

    for (int i = 1; i <= max_iter; i++) {
        TIC;
        if(dcntl[15] != 0) { 
            z = solveM(mu, r);
        } else {
            z = r;
        }
        //IFMASTER clog << " Time in solveM " << TOC << endl;

        TIC;
        rho(0) = dot(r, z);
        
        if (i == 1)
            p = z;
        else {
            beta(0) = rho(0) / rho_1(0);
            p = z + beta(0) * p;
        }
        
        MV_ColMat_double pv(p.ptr(), size_c, 1);
        //IFMASTER clog<< "Time in dot " << TOC << endl;
        TIC;
        Szk = prodSv(pv);
        //IFMASTER clog<< "Time in prodsv " << TOC << endl;
        TIC;

        double *f_ptr = Szk.ptr();
        MV_ColMat_double ff(size_c, 1, 0);
        double *f_o = ff.ptr();

        mpi::all_reduce(inter_comm, f_ptr, size_c, f_o, or_bin);

        q = ff(0);

        alpha(0) = rho(0) / dot(p, q);

        x += alpha(0) * p;
        r -= alpha(0) * q;

        resid = norm(r) / normb;

        if (resid <= tol) {
            tol = resid;
            max_iter = i;
            if(IRANK == 0){
                clog << "Iteration to solve Sz = f " << i << " with a residual of " << resid << endl;
                clog << "Time in iterations : " << MPI_Wtime() - t << endl;
            }
            return x;     
        }
        if(IRANK == 0)
            clog << "Iteration  " << i << " residual " << resid << "\r" << flush;

        //IFMASTER clog<< "Time in otherstuffs " << TOC << endl;
        rho_1(0) = rho(0);
    }
    
    tol = resid;
    if(IRANK == 0){
        clog << "Iteration to solve Sz = f " << max_iter << " with a residual of " << resid << endl;
        clog << "Time in iterations : " << MPI_Wtime() - t << endl;
    }

    return x;
}       /* -----  end of function abcd::pcgS  ----- */
