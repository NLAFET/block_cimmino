// Copyright Institut National Polytechnique de Toulouse (2014) 
// Contributor(s) :
// M. Zenadi <mzenadi@enseeiht.fr>
// D. Ruiz <ruiz@enseeiht.fr>
// R. Guivarch <guivarch@enseeiht.fr>

// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html"

// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 

// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 

// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.

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

#pragma omp parallel for private(rsum)
        for(int r = 0; r < m; ++r) {
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
    #pragma omp parallel for
    for(int k = 0; k < n + 1 ; ++k) {
        a_rp[k]--;
    }
    #pragma omp parallel for
    for(int k = 0; k < nz ; ++k) {
        a_cp[k]--;
    }

    std::vector<double> dc(n, double(1));
    std::vector<double> dr(m, double(1));
    
    // Scale the matrix
    #pragma omp parallel for
    for(int k = 0; k < n; k++) {
        dc[k] = double(1) / dw[k];
        dcol_[k] *= double(1) / dw[k];
    }

    #pragma omp parallel for
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
    int *rp = A.rowptr_ptr();
    int *ci = A.colind_ptr();
    double *v = A.val_ptr();
    
    #pragma omp parallel for
    for ( int i = 0; i < A.dim(0); i++ ) {
        for ( int j = rp[i]; j < rp[i+1]; j++ ) {
            v[j] = drow[i] * v[j]; 
            v[j] = v[j] * dcol[ci[j]];
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
    #pragma omp parallel for
    for ( int i = 0; i < m; i++ ) {
        b(i) = b(i)*drow_[i];
    }

}		/* -----  end of function abcd::scalRhs  ----- */
    void
abcd::diagScaleRhs ( MV_ColMat_double &B)
{
    #pragma omp parallel for
    for ( int i = 0; i < B.dim(0); i++ ) 
        for ( int j = 0; j < B.dim(1); j++ ) 
            B(i,j) = B(i,j)*drow_[i];

}		/* -----  end of function abcd::scalRhs  ----- */
