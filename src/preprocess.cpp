#include "abcd.h"

extern "C"
{
    void mc77id_(int *, double *);
    void mc77ad_(int *job, unsigned *m, unsigned *n, unsigned *nnz, int *jcst, int *irn, double *a,
                 int *iw, int *liw, double *dw, int *ldw, int *icntl, double *cntl,
                 int *info, double *rinfo);
}


int abcd::preprocess()
{
    int ret = 0;

    if (icntl[9] > 0) {
        drow = VectorXd(m);
        dcol = VectorXd(n);

        drow.setOnes();
        dcol.setOnes();
    }


    if (icntl[9] == 2) {
        ///TODO:Use logging
        std::cout << "[-] Scaling with Infinity" << std::endl;

        ret = abcd::scal_matrix(0);

        if (ret != 0)
            return -10;
    }

    abcd::comp_norm();


    if (icntl[9] >= 1) {
        std::cout << "[-] Scaling with Norm 1 & 2" << std::endl;

        ret = abcd::scal_matrix(1);
        ret = abcd::scal_matrix(2);

        if (ret != 0)
            return -11;
    }

    return 0;
}

int abcd::scal_matrix(int norm)
{
    int ldw, liw, nout, job;
    int *iw;
    double *dw;
    int mc77_icntl[10], mc77_info[10];
    double mc77_dcntl[10], mc77_rinfo[10];

    SparseMatrix<double> rmtx(m, m);
    SparseMatrix<double> cmtx(m, m);

    rmtx.reserve(VectorXi::Constant(m, 1));
    cmtx.reserve(VectorXi::Constant(m, 1));

    mc77id_(mc77_icntl, mc77_dcntl);

    if (mc77_icntl[4] == 0)
        ldw = nz;
    else
        ldw = 0;

    if (mc77_icntl[5] == 0) {
        liw = n * 2;
        ldw = ldw + 4 * n;
    } else {
        liw = n;
        ldw = ldw + 2 * n;
    }

    liw = liw * 2;
    ldw = ldw * 2;

    iw = new int[liw];
    dw = new double[ldw];


    // Increment indices to call fortran code
    for (unsigned k = 0; k < n + 1 ; k++) {
        mtx.outerIndexPtr()[k]++;
        mtx.innerIndexPtr()[k]++;
    }
    for (unsigned k = n + 1; k < nz ; k++) {
        mtx.innerIndexPtr()[k]++;
    }
    if (norm < 0) return -12;

    mc77_icntl[9] = 0;
    if (norm == 2) mc77_icntl[6] = 20;

    job = norm;
    mc77ad_(&job, &m, &n, &nz, mtx.outerIndexPtr(), mtx.innerIndexPtr(), mtx.valuePtr(),
            iw, &liw, dw, &ldw, mc77_icntl, mc77_dcntl, mc77_info, mc77_rinfo);

    if (norm == 0)
        cout << "Distance from 1 (norm inf) : " << mc77_rinfo[0] << endl;
    else
        cout << "Distance from 1 (norm " << norm << ") : " << mc77_rinfo[0] << endl;

    // put them back to 0-based for C/C++
    for (unsigned k = 0; k < n + 1 ; k++) {
        mtx.outerIndexPtr()[k]--;
        mtx.innerIndexPtr()[k]--;
    }
    for (unsigned k = n + 1; k < nz ; k++) {
        mtx.innerIndexPtr()[k]--;
    }

    // Scale the matrix
    for (unsigned k = 0; k < m; k++){
        rmtx.insert(k, k) = 1/dw[k];
        // save the scaling meanwhile
        drow(k) *= 1 / dw[k];
    }
    
    for (unsigned k = 0; k < n; k++){
        cmtx.insert(k, k) = 1/dw[k+m];
        // save the scaling meanwhile
        dcol(k) *= 1 / dw[k + m];
    }
    
    mtx = rmtx * mtx * cmtx;

    delete iw, dw;

    return 0;
}

int abcd::comp_norm()
{
    return 0;
}
