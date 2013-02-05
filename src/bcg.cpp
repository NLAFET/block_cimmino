#include <abcd.h>
#include <Eigen/src/misc/blas.h>

void abcd::bcg()
{
    mpi::communicator world;
    // s is the block size of the current run
    int s = std::max<int>(block_size, nrhs);
    if(s < 1) throw - 51;

    //exit(0);
    if(!use_xk) {
        Xk = MV_ColMat_double(n, nrhs, 0);
        //xk = Eigen::MatrixXd(n, nrhs);
        //xk.setZero();
    }

    MV_ColMat_double u(B.dim(0), s, 0);
    MV_ColMat_double p(n, s, 0);
    MV_ColMat_double qp(n, s, 0);
    MV_ColMat_double r(n, s, 0);
    MV_ColMat_double gammak(s, s, 0);
    MV_ColMat_double betak(s, s, 0);
    MV_ColMat_double lambdak(s, nrhs, 0);
    MV_ColMat_double prod_gamma(nrhs, nrhs, 0);
    MV_ColMat_double e1(s, nrhs, 0);


    double *qp_ptr = qp.ptr();
    double *betak_ptr = betak.ptr();
    double *l_ptr = lambdak.ptr();

    double lnrmBs = B.squaredSum();
    double thresh = 1e-12;

    // Sync B norm :
    mpi::all_reduce(inter_comm, &lnrmBs, 1,  &nrmB, std::plus<double>());
    nrmB = sqrt(nrmB);

    for(int k =0; k < e1.dim(1); k++) e1(0,k) = 1;
    char up = 'U';
    char left = 'L';
    char right = 'R';
    char tr = 'T';
    char notr = 'N';
    double alpha = 1;
    bool stay_alive = true;

    // **************************************************
    // ITERATION k = 0                                 *
    // **************************************************

    //for(int k = 0; k < B.dim(1); k++){
        //VECTOR_double vt = B(k);
        //u.setCol(vt, k);
    //}

    mpi::broadcast(intra_comm, stay_alive, 0);

    if(use_xk) {
        MV_ColMat_double sp = sumProject(1e0, B, -1e0, Xk);
        r.setCols(sp, 0, 1);
    } else {
        MV_ColMat_double sp = sumProject(1e0, B, 0, Xk);
        r.setCols(sp, 0, 1);
    }

    
    if(s > nrhs) {
        double *rdata = new double[n * (s - nrhs)];

        srand((unsigned)time(0)); 
        for(int i=0; i< n*(s-nrhs); i++){ 
            rdata[i] = (rand()%n)+1; 
        }

        MV_ColMat_double RR(rdata, n, s-nrhs);
        r.setCols(RR, nrhs, s-nrhs);
    }

    // orthogonalize
    // r = r*gamma^-1
    if(gqr(r, r, gammak, s, false) != 0){
        gmgs(r, r, gammak, s, false);
    }


    p = r;
    prod_gamma = gammak(MV_VecIndex(0, nrhs -1), MV_VecIndex(0, nrhs -1));


    int it = 0;
    double rho = 1;
    std::vector<double> grho(inter_comm.size());
    double mrho;

    double ti = MPI_Wtime();

    if(inter_comm.rank() == inter_comm.size() - 1) {
        cout << "ITERATION " << 0 << endl;
    }
    rho = compute_rho(Xk, B, thresh);

    while((rho > thresh) && (it < itmax)) {
        if(inter_comm.rank() == inter_comm.size() - 1) {
            // a simple hack to clear the screen on unix systems
            cout << "ITERATION " << it + 1 << " Rho = " << rho << endl;
        }
        it++;

        stay_alive = true;
        mpi::broadcast(intra_comm, stay_alive, 0);

        qp = sumProject(0e0, B, 1e0, p);


        if(gqr(p, qp, betak, s, true) != 0){
            gmgs(p, qp, betak, s, true);
        }
        lambdak = e1;

        dtrsm_(&left, &up, &tr, &notr, &s, &nrhs, &alpha, betak_ptr,
               &s, l_ptr, &n);

        lambdak = gemmColMat(lambdak, prod_gamma);


        MV_ColMat_double pl = gemmColMat(p, lambdak);

        Xk(MV_VecIndex(), MV_VecIndex(0, nrhs -1)) += pl(MV_VecIndex(), MV_VecIndex(0, nrhs-1));
        //xk.block(0, 0, n, nrhs) += (p * lambdak).block(0, 0, n, nrhs);
        

        // R = R - QP * B^-T
        dtrsm_(&right, &up, &tr, &notr, &n, &s, &alpha, betak_ptr, &s, qp_ptr, &n);
        r = r - qp;

        double tc = MPI_Wtime();
        if(gqr(r, r, gammak, s, false) != 0)
            gmgs(r, r, gammak, s, false);

        if(inter_comm.rank() == 0)
            cout << "Time : " << MPI_Wtime() - tc << endl << endl;;

        MV_ColMat_double gu = gammak(MV_VecIndex(0, nrhs -1), MV_VecIndex(0, nrhs -1));
        gu = MV_ColMat_double(upperMat(gu));

        prod_gamma = gemmColMat(gu, prod_gamma);
        gammak = upperMat(gammak);
        MV_ColMat_double bu = upperMat(betak);

        //MV_ColMat_double gt = gammak.transpose();
        //cout << bu.dim(1) << " " << gt.dim(1) << endl;

        //betak = gemmColMat(bu, gt);
        //@TODO change to gammak.transpose()
        betak = gemmColMat(bu, gammak);

        p = r + gemmColMat(p, betak);

        rho = abcd::compute_rho(Xk, B, thresh);
        normres.push_back(rho);
        /*
        //mpi::all_gather(inter_comm, rho, grho);
        //mrho = *std::max_element(grho.begin(), grho.end());
        //cout << "ITERATION " << it
                //<< " rho = " << rho << endl;
        */
    }
    stay_alive = false;
    mpi::broadcast(intra_comm, stay_alive, 0);

    if(inter_comm.rank() == 0) {
        //cout << xk << endl;
        cout << "TIME : " << MPI_Wtime() - ti << endl;
    }
}

double abcd::compute_rho(MV_ColMat_double &x, MV_ColMat_double &u, double thresh)
{
    //double nrmX = x.norm();
    int s = x.dim(1);
    MV_ColMat_double R(m, s, 0);
    int pos = 0;
    int ci;
    double gnrmx, nrmXfmX;
    std::vector<double> vnrmx(inter_comm.size());

    for(int p = 0; p < partitions.size(); p++) {
        VECTOR_double compressed_x(partitions[p].dim(1));
        for(int j = 0; j < x.dim(1); j++) {
            compressed_x = VECTOR_double((partitions[p].dim(1)), 0);

            int x_pos = 0;
            for(int i = 0; i < local_column_index[p].size(); i++) {
                int ci = local_column_index[p][i];
                compressed_x(x_pos) = x(ci, j);
                x_pos++;
            }
            //R.col(j).segment(pos, partitions[p].dim(0)) = partitions[p] * compressed_x;
            VECTOR_double vj(R.dim(0));
            vj(MV_VecIndex(pos, pos+partitions[p].dim(0) - 1)) = partitions[p] * compressed_x;
            R.setCol(vj, j);
        }

        pos += partitions[p].dim(0);
    }
    R = u - R;
    //
    double nrmR, nrmX, rho;
    abcd::get_nrmres(x, nrmR, nrmX, nrmXfmX);
    //if(inter_comm.rank() != 0)
        //cout << nrmR << " " << nrmMtx << " " << nrmX << " " << nrmB << endl;
    rho = nrmR / (nrmMtx*nrmX + nrmB);
    //cout << "X -> " << x.col(0).norm() << endl;
    //cout << "X -> " << nrmX << " " << x.col(0).norm() << endl;
    //cout << "R -> " << nrmR << endl;
    //cout << "B -> " << nrmB << endl;
    //cout << "A -> " << nrmA << endl;
    //cout << "M -> " << nrmMtx << endl;
    //return R.col(0).norm() / (nrmA * x.col(0).norm() + u.col(0).norm());
    //cout<< R.col(0).norm() / (nrmA * x.col(0).norm() + nrmB) << endl;
    if(inter_comm.rank() == 0) {
        cout << "Rho = " << rho << endl;
        if(use_xf) cout << "Forward = " << nrmXfmX/nrmXf << endl << endl;
    }
    return rho;
}

void abcd::gmgs(MV_ColMat_double &P, MV_ColMat_double &AP, MV_ColMat_double &R, int s, bool use_a)
{

    int *gr = new int[AP.dim(0)];
    int *gc = new int[AP.dim(0) + 1];
    double *gv = new double[AP.dim(0)];

    for(int i = 0; i < AP.dim(0); i++){
        gr[i] = i;
        gc[i] = i;
        gv[i] = 1;
    }
    gc[AP.dim(0)] = AP.dim(0);

    CompCol_Mat_double G(AP.dim(0), AP.dim(0), AP.dim(0), gv, gr, gc);

    abcd::gmgs(P, AP, R, G, s, use_a);
}

int abcd::gqr(MV_ColMat_double &P, MV_ColMat_double &AP, MV_ColMat_double &R, int s, bool use_a)
{
    int *gr = new int[AP.dim(0)];
    int *gc = new int[AP.dim(0) + 1];
    double *gv = new double[AP.dim(0)];

    for(int i = 0; i < AP.dim(0); i++){
        gr[i] = i;
        gc[i] = i;
        gv[i] = 1;
    }
    gc[AP.dim(0)] = AP.dim(0);

    CompCol_Mat_double G(AP.dim(0), AP.dim(0), AP.dim(0), gv, gr, gc);

    return gqr(P, AP, R, G, s, use_a);
}

/** Generalized modifed Gram-Schmidt
 *
 */
void abcd::gmgs(MV_ColMat_double &p, MV_ColMat_double &ap, MV_ColMat_double &r,
                CompCol_Mat_double g, int s, bool use_a)
{
    r = MV_ColMat_double(s, s, 0);

    // OK!! we have our R here, lets have some fun :)
    for(int k = 0; k < s; k++) {
        VECTOR_double p_k = p(k);
        VECTOR_double ap_k = ap(k);
        r(k, k) = abcd::ddot(p_k, ap_k);

        // if it's negative, make it positive
        if(r(k, k) < r(0,0)*1e-16) {
            r(k, k) = 1;
            cout << "PROBLEM IN GMGS : FOUND AN EXT. SMALL ELMENT " <<  r(k, k) << " postion " << k << endl;
        }
        if(r(k, k) < 0) {
            r(k, k) = abs(r(k, k));
            cout << "PROBLEM IN GMGS : FOUND A NEGATIVE ELMENT " << r(k, k) << " postion " << k << endl;
        }
        if(r(k, k) == 0) {
            r(k, k) = 1;
            cout << "PROBLEM IN GMGS : FOUND A ZERO ELMENT " <<  r(k, k) << " postion " << k << endl;
        }
        // if all is fine, do an sqrt:
        r(k, k) = sqrt(r(k, k));
        p_k = p_k / r(k, k);
        p.setCol(p_k, k);
        if(use_a){
            ap_k = ap_k / r(k, k);
            ap.setCol(ap_k, k);
        }

        if(k < s - 1) {
            for(int j = k + 1; j < s - k + 1; j++){
                ap_k = ap(j);
                p_k = p(k);
                r(k, j) = abcd::ddot(p_k, ap_k);
                VECTOR_double tcol = p(k);
                p_k = p_k - tcol * r(k, j);
                p.setCol(p_k, j);

                if(use_a) {
                    tcol = ap(k);
                    ap_k = ap_k - tcol * r(k, j);
                    ap.setCol(ap_k, j);
                }
            }
            //r.block(k, k + 1, 1, s - (k + 1)) = abcd::ddot(p_k, ap.block(0, k + 1, n, s - (k + 1)));

            //Eigen::MatrixXd tcol = p.col(k);
            //p.block(0, k + 1, n, s - (k + 1))   -=  tcol * r.block(k, k + 1, 1, s - (k + 1));

            //if(use_a) {
                //tcol = ap.col(k);
                //ap.block(0, k + 1, n, s - (k + 1))  -=  tcol * r.block(k, k + 1, 1, s - (k + 1));
            //}
        }
    }
}

int abcd::gqr(MV_ColMat_double &p, MV_ColMat_double &ap, MV_ColMat_double &r,
              CompCol_Mat_double g, int s, bool use_a)
{
    MV_ColMat_double loc_p(n, s, 0);
    MV_ColMat_double loc_ap(n, s, 0);
    MV_ColMat_double loc_r(s, s, 0);

    int pos = 0;
    //  A corriger 
    if(use_a) {
        for(int i = 0; i < n; i++) {
            if(comm_map[i] == 1) {
                for(int j = 0; j < s; j++) {
                    loc_p(pos, j) = p(i, j);
                    loc_ap(pos, j) = ap(i, j);
                }
                pos++;
            }
        }
    } else {
        for(int i = 0; i < p.dim(0); i++) {
            if(comm_map[i] == 1) {
                for(int j = 0; j < p.dim(1); j++) {
                    loc_p(pos, j) = p(i, j);
                }
                pos++;
            }
        }
    }
    int ierr = 0;
    char no = 'N';
    char trans = 'T';
    double alpha, beta;

    alpha = 1;
    beta  = 0;

    double *p_ptr = loc_p.ptr();
    double *ap_ptr = loc_ap.ptr();
    double *l_r_ptr = loc_r.ptr();


    // R = P'AP
    if(use_a){
        int lda_p = loc_p.lda();
        int lda_ap = loc_ap.lda();
        int loc_n = ap.dim(0);

        dgemm_(&trans, &no, &s, &s, &loc_n, &alpha, p_ptr, &lda_p, ap_ptr, &lda_ap, &beta, l_r_ptr, &s);
        //loc_r = MV_ColMat_double(r_ptr, s, s);
    } else{
        int lda_p = loc_p.lda();
        int loc_n = p.dim(0);

        dgemm_(&trans, &no, &s, &s, &loc_n, &alpha, p_ptr, &lda_p, p_ptr, &lda_p, &beta, l_r_ptr, &s);
        //loc_r = MV_ColMat_double(r_ptr, s, s);
    }

    char up = 'U';
    char right = 'R';

    double *r_ptr = r.ptr();
    int rsz = s* s;
    //double tt = MPI_Wtime();
    //mpi::reduce(inter_comm, l_r_ptr, rsz,  r_ptr, std::plus<double>(), 0);
    //mpi::broadcast(inter_comm, r_ptr, rsz , 0);
    mpi::all_reduce(inter_comm, l_r_ptr, s * s,  r_ptr, std::plus<double>());
    //cout << "sub time " << MPI_Wtime() - tt << endl;

    // P = PR^-1    <=> P^T = R^-T P^T
    dpotrf_(&up, &s, r_ptr, &s, &ierr);

//     For the moment if there is an error, just crash!
    if(ierr != 0){
        cout << "PROBLEM IN GQR " << ierr << endl;
        return ierr;
    }
    p_ptr = p.ptr();
    ap_ptr = ap.ptr();

    dtrsm_(&right, &up, &no, &no, &n, &s, &alpha, r_ptr, &s, p_ptr, &n);
    if(use_a){
        dtrsm_(&right, &up, &no, &no, &n, &s, &alpha, r_ptr, &s, ap_ptr, &n);
    }
    //r = MV_ColMat_double(r_ptr, r.dim(0), r.dim(1));

    return 0;
}




