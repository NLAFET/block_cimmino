#include <abcd.h>

extern "C"
{
    void mc77id_(int *, double *);
    void mc77ad_(int *job, int *m, int *n, int *nnz, int *jcst, int *irn, double *a,
                 int *iw, int *liw, double *dw, int *ldw, int *icntl, double *cntl,
                 int *info, double *rinfo);
}

void abcd::scaling()
{
  
    if (m != n) {
        LWARNING << "Matrix is not square, disabling the scaling";
        abcd::icntl[Controls::scaling] = 0;
    }

    if(icntl[Controls::scaling] < 0) {
      dcol_.assign(n, double(1));
      drow_.assign(m, double(1));
    }

    if(icntl[Controls::scaling] >= 0) {
        LINFO << "Row-Scaling with Infinity";

        double rsum;
        dcol_.assign(n, double(1));
        drow_.resize(m);

        for(int r = 0; r < m; r++) {
            rsum = 0;
            for (int c=A.row_ptr(r); c<A.row_ptr(r+1); c++){
                rsum += pow(A(r, A.col_ind(c)), 2);
            }
            drow_[r] = 1/sqrt(rsum);
        }
        abcd::diagScaleMatrix(drow_, dcol_);
    }

    if(icntl[Controls::scaling] >= 1) {
        LINFO << "Scaling with Infinity";
        abcd::scaleMatrix(0);
    }


    if(icntl[Controls::scaling] == 2) {
        LINFO << "Scaling with Norm 1";
        abcd::scaleMatrix(1);

        LINFO << "Scaling with Norm 2";
        abcd::scaleMatrix(2);
    }
}

void abcd::scaleMatrix(int norm)
{
    int ldw, liw, job;
    int *iw;
    double *dw;
    int mc77_icntl[10], mc77_info[10];
    double mc77_dcntl[10], mc77_rinfo[10];

    int *a_rp = A.rowptr_ptr();
    int *a_cp = A.colind_ptr();
    double *a_vp = A.val_ptr();

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
        a_rp[k]++;
    }
    for(int k = 0; k < nz ; k++) {
        a_cp[k]++;
    }

    if(norm < 0){
        delete[] iw;
        delete[] dw;
        throw std::runtime_error("Problem when computing the scaling, got a negative norm.");
    }

    mc77_icntl[6] = 10;

    job = norm;
    mc77ad_(&job, &m, &n, &nz, a_rp, a_cp, a_vp,
            iw, &liw, dw, &ldw, mc77_icntl, mc77_dcntl, mc77_info, mc77_rinfo);


    if(mc77_info[0] < 0) throw - 100 + mc77_info[0];

    if(norm == 0)
        LINFO2 << "Distance from 1 (norm inf) : " << mc77_rinfo[0];
    else
        LINFO2 << "Distance from 1 (norm " << norm << ") : " << mc77_rinfo[0];

    // put them back to 0-based for C/C++
    for(int k = 0; k < n + 1 ; k++) {
        a_rp[k]--;
    }
    for(int k = 0; k < nz ; k++) {
        a_cp[k]--;
    }

    std::vector<double> dc(n, double(1));
    std::vector<double> dr(m, double(1));
    
    // Scale the matrix
    for(int k = 0; k < n; k++) {
        dc[k] = double(1) / dw[k];
        dcol_[k] *= double(1) / dw[k];
    }

    for(int k = 0; k < m; k++) {
        dr[k] = double(1) / dw[k + n];
        drow_[k] *= double(1) / dw[k + n];
    }

    delete[] iw;
    delete[] dw;
    diagScaleMatrix(dr, dc);
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::diagscal
 *  Description:  
 * =====================================================================================
 */
    void
abcd::diagScaleMatrix (std::vector<double> &drow, std::vector<double> &dcol)
{
    for ( int i = 0; i < A.dim(0); i++ ) {
        for ( int j = A.row_ptr(i); j < A.row_ptr(i+1); j++ ) {
            A.val(j) = drow[i] * A.val(j); 
            A.val(j) = A.val(j) * dcol[A.col_ind(j)];
        }
    }
}		/* -----  end of function abcd::diagscal  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::scalRhs
 *  Description:  
 * =====================================================================================
 */
    void
abcd::diagScaleRhs ( VECTOR_double &b)
{
    for ( int i = 0; i < m; i++ ) {
        b(i) = b(i)*drow_[i];
    }

}		/* -----  end of function abcd::scalRhs  ----- */
    void
abcd::diagScaleRhs ( MV_ColMat_double &B)
{
    for ( int i = 0; i < B.dim(0); i++ ) 
        for ( int j = 0; j < B.dim(1); j++ ) 
            B(i,j) = B(i,j)*drow_[i];

}		/* -----  end of function abcd::scalRhs  ----- */
