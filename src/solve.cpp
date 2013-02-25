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
    // if not created yet, do it!
    if( S.dim(0) == 0 ){
        S = abcd::buildS();
    }

    
    /*-----------------------------------------------------------------------------
     *  MUMPS part
     *-----------------------------------------------------------------------------*/
    DMUMPS_STRUC_C mu;
    mu.sym = 0;
    mu.par = 1;
    mu.job = -1;
    mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) inter_comm);

    dmumps_c(&mu);

    mu.icntl[0] = -1;
    mu.icntl[1] = -1;
    mu.icntl[2] = -1;

    if(inter_comm.rank() == 0){ 
        strcpy(mu.write_problem, "/scratch/mzenadi/ss.mtx");
        //mu.icntl[0] = 6;
        //mu.icntl[1] = 6;
        //mu.icntl[2] = 6;
        mu.icntl[3] = 2;
    }

    mu.icntl[7]= 7;

    mu.n = S.dim(0);
    if(inter_comm.size() == 1){ 
        mu.nz= S.NumNonzeros();
        mu.irn= S.rowind_ptr();
        mu.jcn= S.colind_ptr();
        mu.a= S.val_ptr();
        for(int i = 0; i < mu.nz; i++){
            mu.irn[i]++;
            mu.jcn[i]++;
        }
    } else {
        mu.icntl[17]= 3;
        mu.nz_loc = S.NumNonzeros();
        mu.irn_loc = S.rowind_ptr();
        mu.jcn_loc = S.colind_ptr();
        mu.a_loc = S.val_ptr();
        for(int i = 0; i < mu.nz_loc; i++){
            mu.irn_loc[i]++;
            mu.jcn_loc[i]++;
        }
    }
    mu.job = 4;
    dmumps_c(&mu);
    if(mu.info[0] < 0) {
        cout << mu.info[0] << endl;
        exit(0); 
    }
    /*-----------------------------------------------------------------------------
     *  END MUMPS part
     *-----------------------------------------------------------------------------*/

    MV_ColMat_double w;
    if(inter_comm.rank() == 0)
        cout << "[->] Computing w = A^+b" << endl;

    if(dcntl[10] == 0){
        bool stay_alive = true;
        mpi::broadcast(intra_comm, stay_alive, 0);
        w = sumProject(1e0, b, 0e0, Xk);

    } else {
        bcg(b);
        w = Xk; 
    }

    if(inter_comm.rank() == 0)
        cout << "[->] Setting f = - Y w" << endl;

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

    if(inter_comm.rank() == 0){
        mu.rhs = f.ptr();
        mu.nrhs = 1;
    }
    if(inter_comm.rank() == 0)
        cout << "[->] Solving Sz = f" << endl;
    mu.job = 3;
    dmumps_c(&mu);

    if(inter_comm.rank() == 0)
        cout << "[->] Computing zz = (I - P) z" << endl;

    // broadcast f to other cpus, where f is the new z
    double *f_ptr = f.ptr();
    // TODO : better send parts not the whole z
    //
    mpi::broadcast(inter_comm, f_ptr, size_c, 0);

    Xk = MV_ColMat_double(n, 1, 0);
    MV_ColMat_double z_b(m, 1, 0); 
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
        bool stay_alive = true;
        mpi::broadcast(intra_comm, stay_alive, 0);
        f = Xk - sumProject(0e0, b, 1e0, Xk);

    } else {
        //verbose = true;
        use_xk = true;

        if(!use_xk){
            int st = 0;
            for(int p = 0; p < partitions.size(); p++){
                z_b(MV_VecIndex(st, st + partitions[p].dim(0) - 1),
                        MV_VecIndex(0, 0)) = spsmv(partitions[p], local_column_index[p], Xk);
                st += partitions[p].dim(0);
            }

            f = Xk;
        }

        bcg(z_b);

        if(use_xk)
            f = Xk;
        else
            f = f - Xk;
    }


    use_xk = false;

    // the final solution (distributed)
    f = w + f;

    //if(inter_comm.rank()==0){
        //z_b = MV_ColMat_double(m, 1, 0);
        //cout << b.infNorm() << endl;
        //int st = 0;
        //for(int p = 0; p < partitions.size(); p++){
            //z_b(MV_VecIndex(st, st + partitions[p].dim(0) - 1),
                    //MV_VecIndex(0, 0)) = spsmv(partitions[p], local_column_index[p], f);
            //st += partitions[p].dim(0);
        //}
        //z_b = z_b - b;
        //cout << z_b << endl;
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


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::buildS
 *  Description:  
 * =====================================================================================
 */
    Coord_Mat_double
abcd::buildS (  )
{
    Coord_Mat_double shur;
    //MV_ColMat_double S(size_c, size_c, 0);
    std::vector<int> sr, sc;
    std::vector<double> sv;
    bool stay_alive = true;
    if(inter_comm.rank() == 0)
        cout << " [-] Building S" << endl;


    //glob_to_local[n_o + 1]; // just in case it does not exist.

    //std::vector<int> indices;
    //boost::copy(glob_to_local | map_keys, std::back_inserter(indices));

#ifndef NO_USE_XK
    bool use_xk_t = use_xk;
    use_xk = true;
#endif


    std::vector<int> vc, vr;
    std::vector<double> vv;
    //
    //cout << inter_comm.rank() << " getting in" << endl;
    for( int i = 0; i < size_c; i++){
        //if(inter_comm.rank() == 0) cout << "Column : " << i << endl;

        // set xk = [0 I] 
        Xk = MV_ColMat_double(n, 1, 0);
        // Retrieve all keys
        //cout << inter_comm.rank() << " " << indices[0] << " " << indices.back() << " " << n_o + i << endl;
        std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);

        if(iti!=glob_to_local.end()){
            Xk(glob_to_local[n_o + i], 0) = 1;
        }

        MV_ColMat_double b(m, 1, 0); 


        if(dcntl[10] == 0){
            if(iti!=glob_to_local.end()){
                mpi::broadcast(intra_comm, stay_alive, 0);
                MV_ColMat_double sp = Xk - coupleSumProject(0e0, b, 1e0, Xk, size_c);

                //VECTOR_double vv(size_c, 0);
                for(std::map<int,int>::iterator it = glob_to_local.begin(); it != glob_to_local.end(); it++){

                    if(it->first >= n_o){
                        //vv(it->first - n_o) = sp(it->second, 0);
                        vc.push_back(i);
                        vr.push_back(it->first - n_o);
                        vv.push_back(sp(it->second,0));
                    }
                }
                //S.setCol(vv, i);
            } else continue;
        } else {

            //solve with b = 0; and xk_0 = xk
            //if(iti!=glob_to_local.end()){
                //int st = 0;
                //for(int p = 0; p < partitions.size(); p++){
                        //b(MV_VecIndex(st, st + partitions[p].dim(0) - 1),
                                //MV_VecIndex(0, 0)) = spsmv(partitions[p], local_column_index[p], Xk);
                    //st += partitions[p].dim(0);
                //}
            //}
            use_xk = true;
            //verbose = true;
            bcg(b);
        }
        //exit(0);

        //if(inter_comm.rank() == 0){
            //VECTOR_double xx(size_c, 0);
            //for(int k = 0; k < partitions.size(); k++) {
                //for( int j = 0; j < size_c; j++){
                    //if(glob_to_local[n_o+j]!=0){
                        ////cout << n_o + j << " "  << glob_to_local[n_o + j] << " " <<Xk(glob_to_local[n_o + j], 0) << endl;
                    //}
                //}
                ////exit(0);
            //}
        //}
    }
    //cout << inter_comm.rank() << " getting out" << endl;
    //{
        //double *s_ptr = S.ptr();
        //MV_ColMat_double SS(size_c, size_c, 0);
        //double *s_o = SS.ptr();
        //mpi::reduce(inter_comm, s_ptr, size_c*size_c, s_o, std::plus<double>(), 0);
        //S = SS;
    //}
    //
    int *ii = new int[vv.size()];
    int *jj = new int[vv.size()];
    double *val = new double[vv.size()];
    for(int i = 0; i < vv.size(); i++){
        ii[i] = vr[i];
        jj[i] = vc[i];
        val[i] = vv[i];
    }
    shur = Coord_Mat_double(size_c, size_c, vv.size(), val, ii, jj);

    use_xk = false;

    return shur;
}		/* -----  end of function abcd::buildS  ----- */
