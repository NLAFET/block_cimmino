#include "abcd.h"

extern "C"
{
    void mc77id_(int *, double *);
    void mc77ad_(int *job, int *m, int *n, int *nnz, int *jcst, int *irn, double *a,
                 int *iw, int *liw, double *dw, int *ldw, int *icntl, double *cntl,
                 int *info, double *rinfo);
}


void abcd::preprocess()
{
    if(icntl[9] > 0) {
        drow_ = VectorXd(m);
        dcol_ = VectorXd(n);

        drow_.setOnes();
        dcol_.setOnes();
    }

    if(icntl[9] == 1) {
        ///TODO:Use logging
        std::cout << "[-] Scaling with Infinity" << std::endl;

        abcd::scaleMatrix(0);

        SparseMatrix<double> rmtx(m, m);

        rmtx.reserve(VectorXi::Constant(m, 1));
        for(int k = 0; k < m; k++) {
            drow_(k) = sqrt(mtx.row(k).squaredNorm());
            rmtx.insert(k, k) = 1 / drow_[k];
        }
        mtx = rmtx * mtx;
        drow_.setOnes();
        dcol_.setOnes();
        /* END Checking */
    }


    if(icntl[9] == 2) {
        std::cout << "[-] Scaling with Norm 1 & 2" << std::endl;

        abcd::scaleMatrix(1);
        abcd::scaleMatrix(2);
    }

}

void abcd::scaleMatrix(int norm)
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

    if(mc77_icntl[4] == 0)
        ldw = nz;
    else
        ldw = 0;

    if(mc77_icntl[5] == 0) {
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
    for(int k = 0; k < n + 1 ; k++) {
        mtx.outerIndexPtr()[k]++;
        mtx.innerIndexPtr()[k]++;
    }
    for(int k = n + 1; k < nz ; k++) {
        mtx.innerIndexPtr()[k]++;
    }
    if(norm < 0) throw - 12;

    mc77_icntl[9] = 0;
    if(norm == 2) mc77_icntl[6] = 20;

    job = norm;
    mc77ad_(&job, &m, &n, &nz, mtx.outerIndexPtr(), mtx.innerIndexPtr(), mtx.valuePtr(),
            iw, &liw, dw, &ldw, mc77_icntl, mc77_dcntl, mc77_info, mc77_rinfo);

    if(mc77_info[0] < 0) throw - 100 + mc77_info;

    if(norm == 0)
        cout << "Distance from 1 (norm inf) : " << mc77_rinfo[0] << endl;
    else
        cout << "Distance from 1 (norm " << norm << ") : " << mc77_rinfo[0] << endl;

    // put them back to 0-based for C/C++
    for(int k = 0; k < n + 1 ; k++) {
        mtx.outerIndexPtr()[k]--;
        mtx.innerIndexPtr()[k]--;
    }
    for(int k = n + 1; k < nz ; k++) {
        mtx.innerIndexPtr()[k]--;
    }

    // Scale the matrix
    for(int k = 0; k < m; k++) {
        cmtx.insert(k, k) = 1 / dw[k];
        // save the scaling meanwhile
        dcol_(k) *= 1 / dw[k];
    }

    for(int k = 0; k < n; k++) {
        rmtx.insert(k, k) = 1 / dw[k + m];
        // save the scaling meanwhile
        drow_(k) *= 1 / dw[k + m];
    }

    mtx = rmtx * mtx * cmtx;

    delete iw, dw;
}
