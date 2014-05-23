#include <abcd.h>
#include "blas.h"

#include <iostream>
void abcd::bcg(MV_ColMat_double &b)
{
    std::streamsize oldprec = std::cout.precision();

    double t, t1, t2, t1_total, t2_total;
    
    double threshold = dcntl[Controls::threshold];
    int block_size = icntl[Controls::block_size];

    int itmax = icntl[Controls::itmax];
    // s is the block size of the current run
    int s = std::max<int>(block_size, nrhs);

    if (itmax < 0) {
        info[Controls::status] = -11;
        throw std::runtime_error("Max iter number should be at least zero (0)");
    }

    //exit(0);
    if(!use_xk) {
        Xk = MV_ColMat_double(n, nrhs, 0);
    }

    // temporary solution
    MV_ColMat_double u(m, nrhs, 0);
    u = b(MV_VecIndex(0, b.dim(0)-1), MV_VecIndex(0,nrhs-1));

    MV_ColMat_double p(n, s, 0);
    MV_ColMat_double qp(n, s, 0);
    MV_ColMat_double r(n, s, 0);
    MV_ColMat_double gammak(s, s, 0);
    MV_ColMat_double betak(s, s, 0);
    MV_ColMat_double lambdak(s, nrhs, 0);
    MV_ColMat_double prod_gamma(nrhs, nrhs, 0);
    MV_ColMat_double e1(s, nrhs, 0);


    double thresh = threshold;

    double *qp_ptr = qp.ptr();
    double *betak_ptr = betak.ptr();
    double *l_ptr = lambdak.ptr();

    double lnrmBs = infNorm(u);

    // Sync B norm :
    mpi::all_reduce(inter_comm, &lnrmBs, 1,  &nrmB, mpi::maximum<double>());
    mpi::broadcast(inter_comm, nrmMtx, 0);

    for(int k =0; k < e1.dim(1); k++) e1(0,k) = 1;
    char up = 'U';
    char left = 'L';
    char right = 'R';
    char tr = 'T';
    char notr = 'N';
    double alpha = 1;


    // **************************************************
    // ITERATION k = 0                                 *
    // **************************************************

    t1_total = MPI_Wtime();
    if(use_xk) {
        MV_ColMat_double sp = sumProject(1e0, b, -1e0, Xk);
        r.setCols(sp, 0, s);
    } else {
        MV_ColMat_double sp = sumProject(1e0, b, 0, Xk);
        r.setCols(sp, 0, s);
    }


    t1_total = MPI_Wtime() - t1_total;
    // orthogonalize
    // r = r*gamma^-1
    if(gqr(r, r, gammak, s, false) != 0){
        gmgs(r, r, gammak, s, false);
    }

    p = r;
    {
        MV_ColMat_double gu = gammak(MV_VecIndex(0, nrhs -1), MV_VecIndex(0, nrhs -1));
        prod_gamma = MV_ColMat_double(upperMat(gu));
    }


    int it = 0;
    double rho = 1;
    std::vector<double> grho(inter_comm.size());

    double ti = MPI_Wtime();

    t2_total = MPI_Wtime();
    rho = compute_rho(Xk, u);
    t2_total = MPI_Wtime() - t2_total;
    if(comm.rank() == 0) {
        LINFO2 << "ITERATION 0  rho = " << scientific << rho << setprecision(oldprec);
    }

    while(true) {
        it++;
        t = MPI_Wtime();

        qp = sumProject(0e0, b, 1e0, p);

        t1 = MPI_Wtime() - t;

        if(gqr(p, qp, betak, s, true) != 0){
            gmgs(p, qp, betak, s, true);
        }

        lambdak = e1;

        dtrsm_(&left, &up, &tr, &notr, &s, &nrhs, &alpha, betak_ptr, &s, l_ptr, &s);

        lambdak = gemmColMat(lambdak, prod_gamma);


        MV_ColMat_double pl = gemmColMat(p, lambdak);

        Xk(MV_VecIndex(0, Xk.dim(0) - 1), MV_VecIndex(0, nrhs -1)) += 
            pl(MV_VecIndex(0, pl.dim(0)-1), MV_VecIndex(0, nrhs - 1));

        t2 = MPI_Wtime();
        rho = abcd::compute_rho(Xk, u);
        t2 = MPI_Wtime() - t2;
        normres.push_back(rho);

        if((rho < thresh) || (it >= itmax)) break;

        // R = R - QP * B^-T
        dtrsm_(&right, &up, &tr, &notr, &n, &s, &alpha, betak_ptr, &s, qp_ptr, &n);
        r = r - qp;

        if(gqr(r, r, gammak, s, false) != 0)
            gmgs(r, r, gammak, s, false);

        MV_ColMat_double gu = gammak(MV_VecIndex(0, nrhs -1), MV_VecIndex(0, nrhs -1));
        gu = MV_ColMat_double(upperMat(gu));

        prod_gamma = gemmColMat(gu, prod_gamma);
        gammak = upperMat(gammak);
        MV_ColMat_double bu = upperMat(betak);

        // bu * gammak^T
        betak = gemmColMat(bu, gammak, false, true);

        p = r + gemmColMat(p, betak);

        //mpi::all_gather(inter_comm, rho, grho);
        //mrho = *std::max_element(grho.begin(), grho.end());
        //
        t = MPI_Wtime() - t;
        if(comm.rank() == 0 && icntl[Controls::verbose_level] >= 2) {
            int ev = icntl[Controls::verbose_level] >= 3 ? 1 : 10;
            LOG_EVERY_N(ev, INFO) << "ITERATION " << it <<
                " rho = " << scientific << rho <<
                "  Timings: " << setprecision(2) << t <<
                setprecision(oldprec); // put precision back to what it was before
            
        }
        t1_total += t1;
        t2_total += t2;
    }
    

    if(inter_comm.rank() == 0) {
        LINFO2 << "BCG Rho: " << scientific << rho ;
        LINFO2 << "BCG Iterations : " << setprecision(2) << it ;
        LINFO2 << "BCG TIME : " << MPI_Wtime() - ti ;
        LINFO2 << "SumProject time : " << t1_total ;
        LINFO2 << "Rho Computation time : " << t2_total ;
    }
    if (icntl[Controls::aug_type] != 0)
        return;

    info[Controls::nb_iter] = it;

    if(IRANK == 0) {
        sol = MV_ColMat_double(n_o, nrhs, 0);
        solution = sol.ptr();
    }
    
    centralizeVector(solution, nrhs, n_o,
                     Xk.ptr(), nrhs, n);
}

double abcd::compute_rho(MV_ColMat_double &x, MV_ColMat_double &u)
{
    double nrmXfmX;
    double nrmR, nrmX, rho;
    abcd::get_nrmres(x, u, nrmR, nrmX, nrmXfmX);
    rho = nrmR / (nrmMtx*nrmX + nrmB);
    dinfo[Controls::backward] = rho;
    dinfo[Controls::residual] = nrmR;
    dinfo[Controls::scaled_residual] = nrmR/nrmB;
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
    delete[] gr;
    delete[] gc;
    delete[] gv;
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

    delete[] gr;
    delete[] gc;
    delete[] gv;

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

        if(abs(r(k, k)) < abs(r(0,0))*1e-16) {
            r(k, k) = 1;
            LWARNING << "PROBLEM IN GMGS : FOUND AN EXT. SMALL ELMENT " <<  r(k, k) << " postion " << k;
        }
        // if it's negative, make it positive
        if(abs(r(k, k)) < 0) {
            r(k, k) = abs(r(k, k));
            LWARNING << "PROBLEM IN GMGS : FOUND A NEGATIVE ELMENT " << r(k, k) << " postion " << k;
        }
        if(r(k, k) == 0) {
            r(k, k) = 1;
            LWARNING << "PROBLEM IN GMGS : FOUND A ZERO ELMENT " <<  r(k, k) << " postion " << k;
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
            for(int j = k + 1; j < s - k; j++){
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

    } else{
        int lda_p = loc_p.lda();
        int loc_n = p.dim(0);

        dgemm_(&trans, &no, &s, &s, &loc_n, &alpha, p_ptr, &lda_p, p_ptr, &lda_p, &beta, l_r_ptr, &s);
    }



    char up = 'U';
    char right = 'R';

    double *r_ptr = r.ptr();

    mpi::all_reduce(inter_comm, l_r_ptr, s * s,  r_ptr, std::plus<double>());

    // P = PR^-1    <=> P^T = R^-T P^T
    dpotrf_(&up, &s, r_ptr, &s, &ierr);

    if(ierr != 0){
        stringstream err;
        LWARNING << "PROBLEM IN GQR " << ierr << " " << inter_comm.rank();
        LWARNING << "Switching to GMGS";
        
        return ierr;
    }

    p_ptr = p.ptr();
    ap_ptr = ap.ptr();

    dtrsm_(&right, &up, &no, &no, &n, &s, &alpha, r_ptr, &s, p_ptr, &n);
    if(use_a){
        dtrsm_(&right, &up, &no, &no, &n, &s, &alpha, r_ptr, &s, ap_ptr, &n);
    }

    return 0;
}




