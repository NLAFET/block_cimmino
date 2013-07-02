#include <abcd.h>

extern "C"
{
    void mc77id_(int *, double *);
    void mc77ad_(int *job, int *m, int *n, int *nnz, int *jcst, int *irn, double *a,
                 int *iw, int *liw, double *dw, int *ldw, int *icntl, double *cntl,
                 int *info, double *rinfo);
}


void abcd::preprocess()
{
    dcol_ = VECTOR_double(n, double(1));
    drow_ = VECTOR_double(m, double(1));

    if(icntl[9] >= 1) {
        //drow_ = VECTOR_double(m, double(1));

        std::clog << "[-] Scaling with Infinity" << std::endl;

        ///BUG: This takes too much time!
        //
        double rsum;
        drow_ = VECTOR_double(m);
        dcol_ = VECTOR_double(n, double(1));
        for(int r = 0; r < m; r++) {
            rsum = 0;
            for (int c=A.row_ptr(r); c<A.row_ptr(r+1); c++){
                rsum += pow(A(r, A.col_ind(c)), 2);
            }
            drow_(r) = 1/sqrt(rsum);
        }
        abcd::diagScaleMatrix(drow_, dcol_);
        //drow_ = VECTOR_double(m, double(1));
        //dcol_ = VECTOR_double(n, double(1));
        //
        abcd::scaleMatrix(0);

    }


    if(icntl[9] == 2) {
        std::cout << "[-] Scaling with Norm 1" << std::endl;
        abcd::scaleMatrix(1);

        std::cout << "[-] Scaling with Norm 2" << std::endl;
        abcd::scaleMatrix(2);

        //double rsum;
        //VECTOR_double dc_ = VECTOR_double(n, double(1));
        //for(int r = 0; r < m; r++) {
            //rsum = 0;
            //for (int c=A.row_ptr(r); c<A.row_ptr(r+1); c++){
                //rsum += pow(A(r, A.col_ind(c)), 2);
            //}
            //drow_(r) *= 1/sqrt(rsum);
        //}
        //abcd::diagScaleMatrix(drow_, dc_);
        //
        double min_r = 999;
        double max_r = 0;
        double min_c = 999;
        double max_c = 0;

        for(int i = 0; i < m; i++){
            if(min_r > abs(drow_[i])) min_r = abs(drow_[i]);
            if(max_r < abs(drow_[i])) max_r = abs(drow_[i]);
        }
        for(int i = 0; i < n; i++){
            if(min_c > abs(dcol_[i])) min_c = abs(dcol_[i]);
            if(max_c < abs(dcol_[i])) max_c = abs(dcol_[i]);
        }

        cout << "min/max row, col " << min_r << " " << max_r << " | " << min_c << " " << max_c << endl;
    }


    nrmMtx = 0;
    for(int r = 0; r < m; r++) {
        double rsum = 0;
        for (int c=A.row_ptr(r); c<A.row_ptr(r+1); c++){
            rsum += abs(A(r, A.col_ind(c)));
        }
        if(nrmMtx < rsum) nrmMtx = rsum;
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

    if(norm < 0) throw - 12;

    mc77_icntl[6] = 20;

    job = norm;
    mc77ad_(&job, &m, &n, &nz, a_rp, a_cp, a_vp,
            iw, &liw, dw, &ldw, mc77_icntl, mc77_dcntl, mc77_info, mc77_rinfo);


    if(mc77_info[0] < 0) throw - 100 + mc77_info[0];

    if(norm == 0)
        cout << "Distance from 1 (norm inf) : " << mc77_rinfo[0] << endl;
    else
        cout << "Distance from 1 (norm " << norm << ") : " << mc77_rinfo[0] << endl;

    // put them back to 0-based for C/C++
    for(int k = 0; k < n + 1 ; k++) {
        a_rp[k]--;
    }
    for(int k = 0; k < nz ; k++) {
        a_cp[k]--;
    }

    VECTOR_double dc(m, double(1)), dr(n, double(1));
    // Scale the matrix
    for(int k = 0; k < n; k++) {
        dc(k) = double(1) / dw[k];
        dcol_(k) *= double(1) / dw[k];
    }

    for(int k = 0; k < m; k++) {
        dr(k) = double(1) / dw[k + n];
        drow_(k) *= double(1) / dw[k + n];
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
abcd::diagScaleMatrix ( VECTOR_double drow, VECTOR_double dcol)
{
    for ( int i = 0; i < A.dim(0); i++ ) {
        for ( int j = A.row_ptr(i); j < A.row_ptr(i+1); j++ ) {
            A.val(j) = drow(i) * A.val(j); 
            A.val(j) = A.val(j) * dcol(A.col_ind(j));
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
        b(i) = b(i)*drow_(i);
    }

}		/* -----  end of function abcd::scalRhs  ----- */
    void
abcd::diagScaleRhs ( MV_ColMat_double &B)
{
    for ( int i = 0; i < m; i++ ) {
        B(i,0) = B(i,0)*drow_(i);
    }

}		/* -----  end of function abcd::scalRhs  ----- */
