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
        cout << "*      ----------------------      *" << endl;
        cout << "| [--] Building S = Y (I - P) Y^T  |" << endl;
        cout << "*      ----------------------      *" << endl;
    }
    // if not created yet, do it!
    t = MPI_Wtime();
    if( S.dim(0) == 0 ){
        S = abcd::buildS();
    }

    inter_comm.barrier(); //useless! used for timing
    if(inter_comm.rank() == 0){
        cout << "*                                  *" << endl;
        cout << "  [--] T.build S :     " << MPI_Wtime() - t << endl;
        cout << "*                                  *" << endl;
    }

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
        //strcpy(mu.write_problem, "/tmp/sss.mtx");
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

    if(inter_comm.rank() == 0){
        cout << "*                                  *" << endl;
        cout << "  [--] T.Analyse S :   " << MPI_Wtime() - t << endl;
        cout << "*                                  *" << endl;
    }

    t = MPI_Wtime();

    mu.job = 2;
    dmumps_c(&mu);

    if(inter_comm.rank() == 0){
        cout << "*                                  *" << endl;
        cout << "  [--] T.Factorize S : " << MPI_Wtime() - t << endl;
        cout << "*                                  *" << endl;
    }

    if(mu.info[0] < 0) {
        cout << mu.info[0] << endl;
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
        cout << "*      -------------------         *" << endl;
        cout << "| [--] Computing z = S^-1 f        |" << endl;
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
    // broadcast f to other cpus, where f is the new z
    double *f_ptr = f.ptr();
    // TODO : better send parts not the whole z
    //
    mpi::broadcast(inter_comm, f_ptr, size_c, 0);

    return f;
}		/* -----  end of function abcd::solveS  ----- */

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
        for( int i = 0; i < cols.size(); i++){
            std::map<int,int>::iterator iti = glob_to_local.find(n_o + cols[i]);
            if(iti!=glob_to_local.end()) my_cols.push_back(cols[i]);
        }
    }

    if(dcntl[10] == 0 || icntl[15] == 2){

        std::vector<int>::iterator pos = my_cols.begin();
        std::vector<int>::iterator end_pos;

        int share = icntl[14];
        while(pos != my_cols.end()){
            if(pos + share < my_cols.end()) end_pos = pos + share;
            else end_pos = my_cols.end();

            std::vector<int> cur_cols;

            std::copy(pos, end_pos, std::back_inserter(cur_cols));

            int mumps_share = share > 32 ? share/2 : 16;
            setMumpsIcntl(27, mumps_share);
            MV_ColMat_double sp = spSimpleProject(cur_cols);

            for( int j = 0; j < cur_cols.size(); j++){
                int c = cur_cols[j];
                for( int i = 0; i < size_c; i++){
                    if(sp(i,j)!=0){
                        vc.push_back(c);
                        vr.push_back(i);
                        vv.push_back(sp(i, j));
                    }
                }
            }
            pos = end_pos;
        }


    } else {

        //for( int j = 0; j < my_cols.size(); j++){
            //int i = my_cols[j];
        for( int i = 0; i < size_c; i++){
        //
            int my_bro = -1;
            if(inter_comm.rank() == 0) cout << " Column " << i << " out of " << size_c << endl;

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
                //for(int p = 0; p < partitions.size(); p++){
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

                    //cout << "Comm map : " << it->first << " " << comm_map[it->first] << endl;
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

    int *ii = new int[vv.size()];
    int *jj = new int[vv.size()];
    double *val = new double[vv.size()];
    for(int i = 0; i < vv.size(); i++){
        ii[i] = vr[i];
        jj[i] = vc[i];
        val[i] = vv[i];
    }
    shur = Coord_Mat_double(size_c, size_c, vv.size(), val, ii, jj);

    delete[] ii;
    delete[] jj;
    delete[] val;

    use_xk = false;

    return shur;
}		/* -----  end of function abcd::buildS  ----- */



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

    std::map<int,int>::iterator iti;
    for(std::map<int,int>::iterator it = glob_to_local.begin(); it != glob_to_local.end(); it++){

        if(it->first >= n_o){
            for(int j = 0; j < V.dim(1); j++) W(it->second , j) = V(it->first - n_o, j);

        }
    }

    W = W - sumProject(0e0, b, 1e0, W);

    for(std::map<int,int>::iterator it = glob_to_local.begin(); it != glob_to_local.end(); it++){

        if(it->first >= n_o){
            for(int j = 0; j < V.dim(1); j++) R(it->first - n_o, j) = W(it->second , j);

        }
    }

    return R;
}		/* -----  end of function abcd::prodSv  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::buildM
 *  Description:  
 * =====================================================================================
 */
    DMUMPS_STRUC_C 
abcd::buildM (  )
{
    cout << selected_S_columns.size() << endl;
    Coord_Mat_double M = buildS(selected_S_columns);

    S = buildS();

    {
        std::vector<int> mr, mc;
        std::vector<double> mv;

        std::vector<int>::iterator it;
        int r, c;
        double v;
        for(int i = 0; i < M.NumNonzeros(); i++){
            c = M.col_ind(i);
            r = M.row_ind(i);
            v = M.val(i);
            
            if(v == 0) continue;

            if( r != c ) {
                mr.push_back(r);
                mc.push_back(c);
                mv.push_back(v);

                it = find(skipped_S_columns.begin(), skipped_S_columns.end(), r);
                if(it != skipped_S_columns.end()){
                    mr.push_back(c);
                    mc.push_back(r);
                    mv.push_back(v);
                }

            } else if (r == c) {
                mr.push_back(r);
                mc.push_back(c);
                mv.push_back(v);
            }
        }
        cout << "Selected : " << selected_S_columns.size() << endl;
        cout << "Skipped  : " << skipped_S_columns.size() << endl;

        for(int i = 0; i < skipped_S_columns.size(); i++){

            int ro = skipped_S_columns[i];

            mr.push_back(ro);
            mc.push_back(ro);
            mv.push_back(S(ro,ro));

            //double sum = 0;
            //for(int j = 0; j < selected_S_columns.size() ; j++){
                //int co = selected_S_columns[j];
                //sum += M(ro, co) * M(ro, co);
            //}
            //mv.push_back(2*sum);
            //mv.push_back(3);
            //mv.push_back( 0.5*(1 + sqrt(1 - 4*sum)));
            //cout << i << " " << sqrt(1 - 4*sum) << endl;
        }
        
        M = Coord_Mat_double(size_c, size_c, mv.size(), &mv[0], &mr[0], &mc[0]);
    }

    /*-----------------------------------------------------------------------------
     *  MUMPS part
     *-----------------------------------------------------------------------------*/
    double t;
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
        //strcpy(mu.write_problem, "/tmp/mmm.mtx");
        //mu.icntl[0] = 6;
        //mu.icntl[1] = 6;
        //mu.icntl[2] = 6;
        //mu.icntl[3] = 2;
    //}

    mu.n = M.dim(0);

    // parallel analysis if the S is large enough
    //if(mu.n >= 200) {
        //mu.icntl[28 - 1] =  2;
    //}
    mu.icntl[8  - 1] =  7;
    mu.icntl[7  - 1] =  5;
    mu.icntl[14 - 1] =  70;

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

    mu.job = 1;
    dmumps_c(&mu);

    if(inter_comm.rank() == 0){
        cout << "*                                  *" << endl;
        cout << "  [--] T.Analyse M :   " << MPI_Wtime() - t << endl;
        cout << "*                                  *" << endl;
    }

    t = MPI_Wtime();

    mu.job = 2;
    dmumps_c(&mu);

    if(inter_comm.rank() == 0){
        cout << "*                                  *" << endl;
        cout << "  [--] T.Factorize M : " << MPI_Wtime() - t << endl;
        cout << "*                                  *" << endl;
    }

    if(mu.info[0] < 0) {
        cout << mu.info[0] << endl;
        exit(0); 
    }
    /*-----------------------------------------------------------------------------
     *  END MUMPS part
     *-----------------------------------------------------------------------------*/
    return mu;
}		/* -----  end of function abcd::buildM  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::solveM
 *  Description:  
 * =====================================================================================
 */
    VECTOR_double
abcd::solveM (DMUMPS_STRUC_C &mu, VECTOR_double &z )
{
    //VECTOR_double sol(z);
    //if(inter_comm.rank() == 0){
        //mu.rhs = sol.ptr();
        //mu.nrhs = 1;
    //}
    mu.rhs = new double[z.size()];
    for(int i = 0; i < z.size(); i++) mu.rhs[i] = z(i);
    mu.nrhs = 1;

    mu.job = 3;
    dmumps_c(&mu);
    VECTOR_double sol(mu.rhs, z.size());
    delete[] mu.rhs;

    return sol;
}		/* -----  end of function abcd::solveM  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::pcgS
 *  Description:  
 * =====================================================================================
 */
    VECTOR_double
abcd::pcgS ( VECTOR_double &b )
{
    //CG(const Matrix &A, Vector &x, const Vector &b,
   //const Preconditioner &M, int &max_iter, Real &tol)
   //
   DMUMPS_STRUC_C mu = buildM();

   //exit(0);

    double resid, tol = 1e-10;
    int max_iter = size_c;

    VECTOR_double p, z, q;
    VECTOR_double alpha(1), beta(1), rho(1), rho_1(1);

    MV_ColMat_double zk(size_c, 1, 0);
    VECTOR_double x = zk.data();

    double normb = norm(b);
    MV_ColMat_double Szk = prodSv(zk);
    VECTOR_double r = b - Szk(0);

    if (normb == 0.0) 
        normb = 1;
    
    if ((resid = norm(r) / normb) <= tol) {
        tol = resid;
        max_iter = 0;
        return 0;
    }

    for (int i = 1; i <= max_iter; i++) {
        //z = r;
        z = solveM(mu, r);
        rho(0) = dot(r, z);
        //rho(0) = dot(r, r);
        
        if (i == 1)
            p = z;
            //p = r;
        else {
            beta(0) = rho(0) / rho_1(0);
            p = z + beta(0) * p;
            //p = r + beta(0) * p;
        }
        
        //q = A*p;
        MV_ColMat_double pv(p.ptr(), size_c, 1);
        Szk = prodSv(pv);
        q = Szk(0);


        alpha(0) = rho(0) / dot(p, q);
        
        x += alpha(0) * p;
        r -= alpha(0) * q;

        resid = norm(r) / normb;

        if (resid <= tol) {
            tol = resid;
            max_iter = i;
            cout << "Iteration to solve Sz = f " << i << " with a residual of " << resid << endl;
            return x;     
        }

        rho_1(0) = rho(0);
    }
    
    tol = resid;
    cout << "Iteration to solve Sz = f " << max_iter << " with a residual of " << resid << endl;

    return x;
}		/* -----  end of function abcd::pcgS  ----- */
