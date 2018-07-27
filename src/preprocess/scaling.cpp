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
    #ifdef OLD_MUMPS
    void dmumps_simscaleabs_(
            int *irn_loc, int *jcn_loc, double *a_loc,
            int *nz_loc, int *m, int *n, int *numprocs,
            int *myid, int *comm, int *rpartvec, int *cpartvec,
            int *rsndrcvsz, int *csndrcvsz, int *reg,
            int *iwrk, int *iwrksz,
            int *intsz, int *resz, int *op,
            double *rowsca, double *colsca, double *wrkrc, int *iszwrkrc,
            int *sym, int *nb1, int *nb2, int *nb3, double *eps,
            double *onenormerr, double *infnormerr);
    #else
    void dmumps_simscaleabs_(
            int *irn_loc, int *jcn_loc, double *a_loc,
            long long *nz_loc, int *m, int *n, int *numprocs,
            int *myid, int *comm, int *rpartvec, int *cpartvec,
            int *rsndrcvsz, int *csndrcvsz, int *reg,
            int *iwrk, int *iwrksz,
            int *intsz, int *resz, int *op,
            double *rowsca, double *colsca, double *wrkrc, int *iszwrkrc,
            int *sym, int *nb1, int *nb2, int *nb3, double *eps,
            double *onenormerr, double *infnormerr);
    #endif
}

void abcd::scaling()
{
  
    if (m != n) {
        LWARNING << "Matrix is not square, disabling the scaling";
        abcd::icntl[Controls::scaling] = 0;
    }

    dcol_.assign(n, double(1));
    drow_.assign(m, double(1));

    if(icntl[Controls::scaling] !=0) {
        LINFO << "Mumps MC77 scaling with: " << man_scaling[0] << ";" <<
		man_scaling[1] << ";" << man_scaling[2] << ";" << 
		man_scaling[3] << " iterations (#normInf;#norm1;#normInf;#norm2)";

        abcd::scaleMatrix(1);

        diagScaleMatrix(drow_, dcol_);

        if(man_scaling[3] > 0) {
            double rsum;
            std::vector<double> dc(n, double(1));
            std::vector<double> dr(m, double(1));
	    double *Ax = A.val_ptr();

	    #pragma omp parallel for private(rsum)
            for(int r = 0; r < m; ++r) {
                rsum = 0;
                for (int c=A.row_ptr(r); c<A.row_ptr(r+1); c++){
                    rsum += pow(Ax[c], 2);
                }
	        dr[r] = 1/sqrt(rsum);
                drow_[r] *= 1/sqrt(rsum);
            }

            abcd::diagScaleMatrix(dr, dc);
        }
    }
}

void abcd::scaleMatrix(int norm)
{
    int lw, liw, job;
    double eps = 1.e-16;
//    double eps = -1;

    // CSR arrays for the matrix A
    int *a_rp = A.rowptr_ptr();
    int *a_cp = A.colind_ptr();
    double *a_vp = A.val_ptr();
    // Increment column indices to call fortran code
    #pragma omp parallel for
        for(int k = 0; k < nz ; k++) {
            a_cp[k]++;
        }

    // Local scaling arrays
    std::vector<double> dc(n, double(1));
    std::vector<double> dr(m, double(1));

    int issym;
    int numprocs = 1, myid = 0;
    std::vector<int> iwk;
    std::vector<double> dwk;

    int maxmn, intsz, resz;
    int reg[12];
    int *rpartvec, *cpartvec, *rsndrcvsz, *csndrcvsz;

    // Convert row indices from CSR to coordinate format
    // We do not need to increment the array as it won't be used used as is
    std::vector<int> rp(nz);
    for (int i = 0; i < m; ++i) {
        for (int j = a_rp[i]; j < a_rp[i+1]; ++j) {
            rp[j] = i+1;
        }
    }

    maxmn = n > m ? n : m;
    liw = 4 * maxmn;
    iwk.assign(liw, 0);
    rpartvec = new int[m];
    cpartvec = new int[n];

    // hard coded 2*numprocs, we run it in sequential
    rsndrcvsz = new int[2*numprocs];
    csndrcvsz = new int[2*numprocs];

    int co = MPI_Comm_c2f((MPI_Comm) comm);

    // NB1, NB2, NB3: algo runs successively
    // NB1 iters of inf-norm
    // NB2 iters of 1-norm
    // NB3 iters of inf-norm
    int nb1 = man_scaling[0], nb2 = man_scaling[1], nb3 = man_scaling[2];
    issym = 0;
    lw = 0;
    double err, errinf;

    // 1. estimate memory 
    job = 1;
    #ifdef OLD_MUMPS
        int nz_tmp=nz;
    #else
        long long nz_tmp=nz;
    #endif

    dmumps_simscaleabs_(
        &rp[0], a_cp, a_vp, &nz_tmp, &m, &n, &numprocs, &myid, &co,
        rpartvec, cpartvec, rsndrcvsz, csndrcvsz, reg,
        &iwk[0], &liw,
        &intsz, &resz, &job,
        &dr[0], &dc[0], &dwk[0], &lw,
        &issym, &nb1, &nb2, &nb3, &eps, &err, &errinf);

    if (liw < intsz) {
        liw = intsz;
        iwk.assign(liw, 0);
    }

    lw = resz;
    dwk.assign(lw, 0);

    // 2. compute drow and dcol
    job = 2;
    dmumps_simscaleabs_(
        &rp[0], a_cp, a_vp, &nz_tmp, &m, &n, &numprocs, &myid, &co,
        rpartvec, cpartvec, rsndrcvsz, csndrcvsz, reg,
        &iwk[0], &liw,
        &intsz, &resz, &job,
        &dr[0], &dc[0], &dwk[0], &lw,
        &issym, &nb1, &nb2, &nb3, &eps, &err, &errinf);

    // Set actual scaling arrays (directly the inverse)
    for(int k = 0; k < n; k++) {
        dcol_[k] *= dc[k];
    }
    for(int k = 0; k < m; k++) {
        drow_[k] *= dr[k];
    }

    // Free memory
    delete[] rpartvec;
    delete[] cpartvec;
    delete[] rsndrcvsz;
    delete[] csndrcvsz;
 
    #pragma omp parallel for
    for(int k = 0; k < nz ; k++) {
        a_cp[k]--;
    }
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
