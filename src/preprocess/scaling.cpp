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
}

void abcd::scaling()
{
  
    if (m != n) {
        LWARNING << "Matrix is not square, disabling the scaling";
        abcd::icntl[Controls::scaling] = 0;
    }

    dcol_.assign(n, double(1));
    drow_.assign(m, double(1));

    if(icntl[Controls::scaling] > 0) {
        LINFO << "Scaling the matrix";

        abcd::scaleMatrix(0);

        // disable col scale temp.
        dcol_.assign(n, double(1));

        diagScaleMatrix(drow_, dcol_);
    }
}

void abcd::scaleMatrix(int norm)
{
    int lwk, liwk, issym;
    int numprocs = 1, myid = 0, job;
    std::vector<int> iwk;
    std::vector<double> dwk;

    int maxmn, intsz, resz;
    int reg[12];
    int *rpartvec, *cpartvec, *rsndrcvsz, *csndrcvsz;

    int *a_rp = A.rowptr_ptr();
    std::vector<int> rp(nz);
    for (int i = 0; i < m; ++i) {
        for (int j = a_rp[i]; j < a_rp[i+1]; ++j) {
            rp[j] = i;
        }
    }

    int *a_cp = A.colind_ptr();
    double *a_vp = A.val_ptr();

    maxmn = n > m ? n : m;

    liwk = 4 * maxmn;

    iwk.assign(liwk, 0);
    rpartvec = new int[m];
    cpartvec = new int[n];

    // hard coded 2*numprocs, we run it in sequential
    rsndrcvsz = new int[2];
    csndrcvsz = new int[2];

    int co = MPI_Comm_c2f((MPI_Comm) comm);

    // NB1, NB2, NB3: algo runs successively
    // NB1 iters of inf-norm 
    // NB2 iters of 1-norm   
    // NB3 iters of inf-norm 
    int nb1 = 5, nb2 = 10, nb3 = 5;
    double eps = 1.e-8;
    issym = 0;
    lwk = 0;
    double err, errinf;

    std::vector<double> dc(n, 1);
    std::vector<double> dr(m, 1);

    // estimate memory 
    job = 1;
    dmumps_simscaleabs_(
        &rp[0], a_cp, a_vp, &nz, &m, &n, &numprocs, &myid, &co,
        rpartvec, cpartvec, rsndrcvsz, csndrcvsz, reg,
        &iwk[0], &liwk,
        &intsz, &resz, &job,
        &dr[0], &dc[0], &dwk[0], &lwk,
        &issym, &nb1, &nb2, &nb3, &eps, &err, &errinf);

    if (liwk < intsz) {
        liwk = intsz;
        iwk.assign(liwk, 0);
    }

    lwk = resz;
    dwk.assign(lwk, 0);

    // compute drow and dcol
    job = 2;
    dmumps_simscaleabs_(
        &rp[0], a_cp, a_vp, &nz, &m, &n, &numprocs, &myid, &co,
        rpartvec, cpartvec, rsndrcvsz, csndrcvsz, reg,
        &iwk[0], &liwk,
        &intsz, &resz, &job,
        &dr[0], &dc[0], &dwk[0], &lwk,
        &issym, &nb1, &nb2, &nb3, &eps, &err, &errinf);

    for(int k = 0; k < n; k++) {
        dcol_[k] = dc[k];
    }

    for(int k = 0; k < m; k++) {
        drow_[k] = dr[k];
    }

    delete[] rpartvec;
    delete[] cpartvec;
    delete[] rsndrcvsz;
    delete[] csndrcvsz;
    
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
