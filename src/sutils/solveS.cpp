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
    if( S.dim(0) == 0 ){
        S = abcd::buildS();
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
    //
    /*
    mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) inter_comm);
    */

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
        //mu.icntl[3] = 2;
    }

    mu.n = size_c;

    // parallel analysis if the S is large enough
    //if(mu.n >= 200) {
        //mu.icntl[28 - 1] =  2;
    //}
    mu.icntl[8  - 1] =  77;
    mu.icntl[7  - 1] =  5;
    //mu.icntl[6  - 1] =  5;
    //mu.icntl[12 - 1] =  2;
    mu.icntl[14 - 1] =  70;
    mu.icntl[40 - 1] =  1; // no type 2 if 1


    /*
    mu.keep[24 -1 ] = 8;
    mu.keep[78 -1 ] = 0;
    mu.keep[77 -1 ] = 0;
    mu.keep[68 -1 ] = 0;
    mu.keep[9 -1 ] = 4000;
    */

#ifdef CENTRALIZE
    /*  This part is in case it's centralization */

    //get the total nz
    int loc_nz = S.NumNonzeros();
    mpi::all_reduce(inter_comm, loc_nz, mu.nz, std::plus<int>());

    // get everything to the master
    std::vector<int> irn_loc (loc_nz);
    std::vector<int> jcn_loc (loc_nz);
    std::vector<double> a_loc(loc_nz);

    for(int i = 0; i < loc_nz; i++){
#ifdef CENTRALIZED_SUM
        irn_loc[i] = S.row_ind(i);
        jcn_loc[i] = S.col_ind(i);
#else
        irn_loc[i] = S.row_ind(i) + 1;
        jcn_loc[i] = S.col_ind(i) + 1;
#endif
        a_loc[i] = S.val(i);
    }
    
    if(IRANK == 0){

        mu.irn = new int[mu.nz];
        mu.jcn = new int[mu.nz];
        mu.a   = new double[mu.nz];

        memcpy( mu.irn, &irn_loc[0], loc_nz * sizeof(int));
        memcpy( mu.jcn, &jcn_loc[0], loc_nz * sizeof(int));
        memcpy( mu.a  , &a_loc[0],   loc_nz * sizeof(double));

        int current_pos = loc_nz;

        for(int i = 1; i < inter_comm.size(); i++){

            int rnz; 
            inter_comm.recv(i, 70, rnz);

            inter_comm.recv(i, 71, mu.irn + current_pos, rnz);
            inter_comm.recv(i, 72, mu.jcn + current_pos, rnz);
            inter_comm.recv(i, 73, mu.a   + current_pos, rnz);

            current_pos += rnz;
        }

#ifdef CENTRALIZED_SUM
        Coord_Mat_double SS(size_c, size_c, mu.nz, mu.a, mu.irn, mu.jcn);
        CompCol_Mat_double SC(SS);

        
        delete[] mu.irn, mu.jcn, mu.a;
        std::vector<int> ii, jj;
        std::vector<double> vv;

        int r;
        for(int c = 0; c < SC.dim(1); c++){
            std::map<int,double> mi;
            for(int i=SC.col_ptr(c); i < SC.col_ptr(c+1); i++){
                r = SC.row_ind(i);
                if(mi[r]){
                    mi[r] += SC.val(i);
                } else {
                    mi[r] = SC.val(i);
                }
            }
            for(std::map<int,double>::iterator it = mi.begin();
                    it != mi.end(); it++) {
                ii.push_back(it->first);
                jj.push_back(c);
                vv.push_back(it->second);
            }
        }
        mu.nz  = vv.size();
        mu.irn = new int[mu.nz];
        mu.jcn = new int[mu.nz];
        mu.a   = new double[mu.nz];

        for(int i = 0; i < mu.nz; i++){
            mu.irn[i] = ii[i] + 1;
            mu.jcn[i] = jj[i] + 1;
            mu.a  [i] = vv[i] + 1;
        }
#endif

    } else {
        inter_comm.send(0, 70, loc_nz);
        inter_comm.send(0, 71, &irn_loc[0], loc_nz);
        inter_comm.send(0, 72, &jcn_loc[0], loc_nz);
        inter_comm.send(0, 73, &a_loc[0], loc_nz);
    }

    //if(IRANK == 0) clog << "TOTAL " << mu.nz << endl;
    //IBARRIER;
    //exit(0);

#else
    if(inter_comm.size() == 1){ 
        mu.nz= S.NumNonzeros();
        mu.irn = new int[mu.nz];
        mu.jcn = new int[mu.nz];
        mu.a= S.val_ptr();

        for(int i = 0; i < mu.nz; i++){
            mu.irn[i] = S.row_ind(i) + 1;
            mu.jcn[i] = S.col_ind(i) + 1;
        }

    } else {
        mu.icntl[18 - 1]= 3;
        mu.nz_loc = S.NumNonzeros();

        mu.irn_loc = new int[mu.nz_loc];
        mu.jcn_loc = new int[mu.nz_loc];
        mu.a_loc   = new double[mu.nz_loc];
        //mu.irn_loc = S.rowind_ptr();
        //mu.jcn_loc = S.colind_ptr();
        //mu.a_loc = S.val_ptr();
        for(int i = 0; i < mu.nz_loc; i++){
            //mu.irn_loc[i]++;
            //mu.jcn_loc[i]++;
            mu.irn_loc[i] = S.row_ind(i) + 1;
            mu.jcn_loc[i] = S.col_ind(i) + 1;
            mu.a_loc[i] = S.val(i);
        }
    }
#endif

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
