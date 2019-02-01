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

/*!
 * \file splib_utils.cpp
 * \brief Some additions to SparseLib++
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 * \date 01/23/2013 06:47:59 PM
 * Revision:  none
 * Compiler:  gcc
 */

#include <abcd.h>
#include "blas.h"

#define MKL_INT int

//using namespace boost::lambda;

/*!
 *  \brief Compute infinite norm of a vector
 *
 *  Compute infinite norm of a vector
 *
 *  \param V: the vector which norm must be computed
 *
 *  \return the infinite-norm of the vector
 *
 */
double infNorm(VECTOR_double &V){
    double max = 0;
    for (int i = 0; i < V.size(); i++){
        //if(abs(V(i)) >= max) max = abs(V(i));
        max = abs(V(i)) > max ? abs(V(i)) : max;
    }
    return max;
}               /* -----  end of function infNorm  ----- */

/*!
 *  \brief Compute infinite norm of a matrix
 *
 *  Compute infinite norm of a matrix
 *
 *  \param V: the matrix which norm must be computed
 *
 *  \return the infinite-norm of the matrix
 *
 */
double infNorm(MV_ColMat_double &V){
    double max = 0;
    double *v_ptr = V.ptr();
    for (int i = 0; i < V.dim(0)*V.dim(1); i++){
        if(abs(v_ptr[i]) >= max) max = abs(v_ptr[i]);
    }
    return max;
}               /* -----  end of function infNorm  ----- */

/*!
 *  \brief Compute square of the 2-norm of a subvector
 *
 *  Compute square of the 2-norm of a subvector
 *
 *  \param V: the vector which norm of a subpart must be computed
 *  \param I: indices of the elements of interest of the vector
 *
 *  \return the squared 2-norm
 *
 */
double squaredNorm(VECTOR_double &V, VECTOR_int &I){
    double sum = 0;
    double e;
    for (int i = 0; i < I.size(); i++){
        e = V(I(i));
        sum += e*e;
    }
    return sum;
}               /* -----  end of function squaredNorm  ----- */

CompRow_Mat_double CSC_extractByIndices(CompRow_Mat_double &M, vector<int> row_index)
{
  int *mptr = M.rowptr_ptr();
  int * sub_row_vect = new int[row_index.size()+1];
  sub_row_vect[0]=0;
  for(int i = 1; i < row_index.size()+1; i++){
    sub_row_vect[i] = sub_row_vect[i-1] +  mptr[ row_index[i-1]+1 ] -  mptr[ row_index[i-1] ];
    //cout << " subrow vect " << i << " "  << sub_row_vect[i] << " " << row_index[i] << " " << mptr[ row_index[i] ] << " " << mptr[ row_index[i-1]] << endl;
  }
  int local_nnz = sub_row_vect[ row_index.size()];

  int * sub_col_vect =  new int[local_nnz];
  double * sub_val_vect = new double[local_nnz] ;

  int * colptr = M.colind_ptr();
  double * valptr = M.val_ptr();

  int cnt =0;
  for(int i = 0; i < row_index.size(); i++){
    for(int j= mptr[ row_index[i] ]; j < mptr[ row_index[i]+1 ] ; j++){
      sub_col_vect[cnt] = colptr[j];
      sub_val_vect[cnt] =  valptr[j];
      cnt++;
    }
  }

  CompRow_Mat_double nM( row_index.size(), M.dim(1), local_nnz,
      sub_val_vect,
      sub_row_vect,
      sub_col_vect
      );

  //delete[] sub_col_vect, sub_val_vect, m_row_ptr, sub_row_vect;
  std::cout << "Returned submatrix" << std::endl;
  std::cout << nM;
  return nM;
}

// Assume col_index array is sorted and contiguous
CompRow_Mat_double CSR_extractByColIndices(CompRow_Mat_double &M, vector<int> col_index)
{
  int *mptr       = M.rowptr_ptr();
  int * colptr    = M.colind_ptr();
  double * valptr = M.val_ptr();

  int nrow = M.dim(0);
  int col_shift = col_index[0];

  std::cout << M.rowptr_.size() << std::endl;
  printf("Allocate array of size %d\n", nrow + 1);
  int * sub_row_vect = new int[nrow + 1];

  std::cout << "Col_INDEX" <<std::endl;
  for (auto elt = col_index.begin(); elt != col_index.end(); elt++)
    std::cout << *elt <<std::endl;

  sub_row_vect[0] = 0;
  int nlcol = 0; //#column for the current row
  int k     = 0; //pos in col_index
  for(int i = 0; i < nrow; i++){
    nlcol = 0;
    k     = 0;
  //for(int j = mptr[i]; j < mptr[i + 1] ; j++){
    for(int j = mptr[i]; j < mptr[i + 1] ; j++){
      while(k < col_index.size()){
        if (colptr[j] > col_index[k]){
          k++;
        }else{
          if(colptr[j] == col_index[k]){
            printf("row[%d] col %d ? %d => %d\n", i, colptr[j], col_index[k], nlcol+1);
            nlcol++;
          }
          break;
        }
      }
    }
    printf("Want to set %d using %d and %d\n", i+1, i, nlcol);
    sub_row_vect[i + 1] = sub_row_vect[i] + nlcol;
  }
  int local_nnz = sub_row_vect[nrow];

  int * sub_col_vect    = new int[local_nnz];
  double * sub_val_vect = new double[local_nnz];

  int cnt =0;
  for(int i = 0; i < nrow; i++){
    k = 0;
    for(int j = mptr[i]; j < mptr[i + 1] ; j++){
      while(k < col_index.size()){
        if (colptr[j] > col_index[k]){
          k++;
        }else{
          if(colptr[j] == col_index[k]){
            sub_col_vect[cnt] = colptr[j] - col_shift; // /!\  shifted assuming contiguous columns indices in col_index 
            assert(colptr[j] - col_shift >= 0);
            sub_val_vect[cnt] = valptr[j];
            cnt++;
          }
          break;
        }
      }
    }
  }
  assert(cnt == local_nnz);

  for (int i = 0; i < nrow + 1; i++)
    printf("%d ", sub_row_vect[i]);
  printf("\n");
  for (int i = 0; i < local_nnz; i++)
    printf("%d ", sub_col_vect[i]);
  printf("\n");

  printf("Create CSR matrix of size %dx%ld\n", nrow, col_index.size());
  CompRow_Mat_double nM(nrow, col_index.size(), local_nnz,
      sub_val_vect,
      sub_row_vect,
      sub_col_vect
      );

  //delete[] sub_col_vect, sub_val_vect, m_row_ptr, sub_row_vect;
  std::cout << "Returned submatrix" << std::endl;
  std::cout << nM;
  return nM;
}

/*!
 *  \brief Extract rows from a CSR row matrix
 *
 *  Extract rows from a CSR row matrix
 *
 *  \param M: matrix which subpart must be extracted
 *  \param st_row: starting row of the subpart to extract
 *  \param nb_rows: number of rows to extract
 *
 *  \return the extracted matrix
 *
 */
CompRow_Mat_double CSR_middleRows (CompRow_Mat_double &M, int st_row, int nb_rows) {
    int st_index, ed_index;

    // starting index in JCN
    st_index = M.row_ptr(st_row);
    // last index in JCN
    ed_index = M.row_ptr(st_row + nb_rows) - 1;

    int starting_point = M.row_ptr(st_row);

    int * m_row_ptr = M.rowptr_ptr() + st_row;
    int * sub_row_vect = new int[nb_rows + 1];

    for(int i = 0; i <= nb_rows; i++){
        sub_row_vect[i] = m_row_ptr[i] - starting_point;
    }

    int * sub_col_vect = M.colind_ptr() + st_index;
    double * sub_val_vect = M.val_ptr() + st_index;

    CompRow_Mat_double nM( nb_rows, M.dim(1), ed_index - st_index + 1,
            sub_val_vect,
            sub_row_vect,
            sub_col_vect
            );

    //delete[] sub_col_vect, sub_val_vect, m_row_ptr, sub_row_vect;
    return nM;
}               /* -----  end of function CSR_middleRows  ----- */

/*!
 *  \brief Extract rows from a CSC row matrix
 *
 *  Extract rows from a CSC row matrix
 *
 *  \param M: matrix which subpart must be extracted
 *  \param st_row: starting row of the subpart to extract
 *  \param nb_rows: number of rows to extract
 *
 *  \return the extracted matrix
 *
 */
CompCol_Mat_double CSC_middleRows (CompRow_Mat_double &M, int st_row, int nb_rows) {
    return CompCol_Mat_double(CSR_middleRows(M, st_row, nb_rows));
}               /* -----  end of function CSC_middleRows  ----- */

/*!
 *  \brief Extract a column from a matrix
 *
 *  Extract column from a matrix
 *
 *  \param M: matrix which subpart must be extracted
 *  \param col_num: index of the column to extract
 *  \param ind: index of the rows in the column
 *
 *  \return the extracted vector column
 *
 */
VECTOR_double middleCol(CompCol_Mat_double &M, int col_num, VECTOR_int &ind){
    int st_index, ed_index;

    VECTOR_double c(M.dim(0), 0);

    // starting index in JCN
    st_index = M.col_ptr(col_num);
    // last index in JCN
    ed_index = M.col_ptr(col_num + 1) - 1;


    if(st_index > ed_index) return c;

    VECTOR_int sub_row_vect(ed_index - st_index + 1);
    sub_row_vect() = M.row_ind(MV_VecIndex(st_index, ed_index));
    ind = sub_row_vect();

    VECTOR_double sub_val_vect(ed_index - st_index + 1);
    sub_val_vect() = M.val(MV_VecIndex(st_index, ed_index));


    for(int i = 0; i < sub_val_vect.size(); i++){
        c(sub_row_vect(i)) = sub_val_vect(i);
    }

    return c;
}               /* -----  end of function middleCol  ----- */

/*!
 *  \brief Extract columns from a matrix
 *
 *  Extract columns from a matrix
 *
 *  \param M: matrix which subpart must be extracted
 *  \param ci: indices of the column to extract
 *
 *  \return the extracted submatrix
 *
 */
    CompCol_Mat_double
sub_matrix ( CompCol_Mat_double &M, std::vector<int> &ci )
{
    CompCol_Mat_double SM;
    int st_col, ed_col;
    int c = 0, j = 0, nzc;

    int *M_col_ptr = M.colptr_ptr();
    int *M_row_ptr = M.rowind_ptr();
    double *M_val_ptr = M.val_ptr();


    // Compute the number of non-zeros
    for(int k=0; k < ci.size(); ++k){
        st_col = M_col_ptr[ci[k]];
        ed_col = M_col_ptr[ci[k] + 1];

        nzc = ed_col - st_col;

        // if the column is empty
        if (st_col == ed_col) continue;

        c += nzc;
        ++j;
    }

    std::vector<int> v_sm_c(ci.size() + 1);
    std::vector<int> v_sm_r(c);
    std::vector<double> v_sm_v(c);

    c = 0; j = 0;

    // for all columns in ci
    for(int k=0; k < ci.size(); k++){
        st_col = M.col_ptr(ci[k]);
        ed_col = M.col_ptr(ci[k] + 1);

        nzc = ed_col - st_col;

        // if the column is empty
        if (st_col == ed_col) continue;

        v_sm_c[j] = c;
        //
        int pos = st_col;
        for(int i = c; i < c + nzc; ++i){
            v_sm_r[i] = M_row_ptr[pos];
            v_sm_v[i] = M_val_ptr[pos];
            pos++;
        }

        c += nzc;
        j++;
    }
    if (c==0) return SM;
    v_sm_c[j] = c;
    SM = CompCol_Mat_double(M.dim(0), j, c,
            &v_sm_v[0], &v_sm_r[0], &v_sm_c[0]);

    return SM;
}		/* -----  end of function sub_matrix  ----- */


/*!
 *  \brief Compute sparse Matrix-matrix product for row matrices
 *
 *  Compute sparse Matrix-matrix product for row matrices
 *
 *  \param A: first sparse matrix
 *  \param B: second sparse matrix
 *
 *  \return the product matrix
 *
 */
    CompRow_Mat_double
spmm ( CompRow_Mat_double &A, CompRow_Mat_double &BT )
{
    //CompRow_Mat_double BT = csr_transpose(B);

    std::map<int,int> iw;
    std::map<int,int> jc, ic;
    std::map<int,double> c;
    ic[0] = 0;

    double scal;
    int i, ka, kb, j, jcol, len=0;

    for( i = 0; i < A.dim(0); i++){
        for ( ka = A.row_ptr(i); ka < A.row_ptr(i+1); ka++){
            scal = A.val(ka);
            j = A.col_ind(ka);

            for(kb = BT.row_ptr(j); kb < BT.row_ptr(j+1); kb++){
                // i-th row from A x jcol-th column from B
                jcol = BT.col_ind(kb);
                if(iw.find(jcol)==iw.end()){
                    // first time!
                    jc[len] = jcol;
                    c[len]  = scal*BT.val(kb);
                    iw[jcol] = len;
                    len++;
                } else {
                    c[iw[jcol]] += scal*BT.val(kb);

                    if(c[iw[jcol]] == 0){
                        c.erase(iw[jcol]);
                        jc.erase(iw[jcol]);
                        iw.erase(jcol);
                        len--;
                    }
                }
            }
        }
        iw.clear();
        ic[i+1] = len;
    }
    VECTOR_int vc(len), vr(A.dim(0) + 1);
    VECTOR_double vv(len);

    for(i = 0; i < len; i++){
        vc[i] = jc[i];
        vv[i] = c[i];
    }

    for(i = 0; i < A.dim(0) + 1; i++){
        vr[i] = ic[i];
    }

    CompRow_Mat_double C(A.dim(0), BT.dim(1), len, vv, vr, vc);

    return C;
}		/* -----  end of function spmm_original  ----- */


/*!
 *  \brief Compute the transpose of a CSC column matrix
 *
 *  Compute the transpose of a CSC column matrix
 *
 *  \param M: matrix to transpose
 *
 *  \return the transposed matrix
 *
 */
    CompCol_Mat_double
csc_transpose ( CompCol_Mat_double &M )
{
    CompRow_Mat_double Mt = csr_transpose(M);
    return CompCol_Mat_double(Mt);

}		/* -----  end of function csc_transpose  ----- */

/*!
 *  \brief Compute the transpose of a CSC row matrix
 *
 *  Compute the transpose of a CSC row matrix
 *
 *  \param M: matrix to transpose
 *
 *  \return the transposed matrix
 *
 */
    CompCol_Mat_double
csc_transpose ( CompRow_Mat_double &M )
{
    VECTOR_int vc(M.dim(0) + 1);
    VECTOR_int vr(M.NumNonzeros());
    VECTOR_double vv(M.NumNonzeros());

    vc() = M.row_ptr(MV_VecIndex());
    vr() = M.col_ind(MV_VecIndex());
    vv() = M.val(MV_VecIndex());


    return CompCol_Mat_double(M.dim(1), M.dim(0), M.NumNonzeros(), vv, vr, vc);
}		/* -----  end of function csr_transpose  ----- */

/*!
 *  \brief Compute the transpose of a CSR column matrix
 *
 *  Compute the transpose of a CSR column matrix
 *
 *  \param M: matrix to transpose
 *
 *  \return the transposed matrix
 *
 */
    CompRow_Mat_double
csr_transpose ( CompCol_Mat_double &M )
{
    VECTOR_int vr(M.dim(1) + 1);
    VECTOR_int vc(M.NumNonzeros());
    VECTOR_double vv(M.NumNonzeros());

    vr() = M.col_ptr(MV_VecIndex());
    vc() = M.row_ind(MV_VecIndex());
    vv() = M.val(MV_VecIndex());


    return CompRow_Mat_double(M.dim(1), M.dim(0), M.NumNonzeros(), vv, vr, vc);
}		/* -----  end of function csr_transpose  ----- */

/*!
 *  \brief Compute the transpose of a CSR row matrix
 *
 *  Compute the transpose of a CSR row matrix
 *
 *  \param M: matrix to transpose
 *
 *  \return the transposed matrix
 *
 */
    CompRow_Mat_double
csr_transpose ( CompRow_Mat_double &M )
{
    return CompRow_Mat_double(csc_transpose(M));
}		/* -----  end of function csr_transpose  ----- */


/*!
 *  \brief Resize matrix by adding empty columns at the end
 *
 *  Resize matrix by adding empty columns at the end
 *
 *  \param M: matrix to resize
 *  \param new_size: new size of the matrix
 *
 *  \return the resized matrix
 *
 */
    CompCol_Mat_double
resize_columns ( CompCol_Mat_double &M, int new_size )
{
    if(new_size < M.dim(1)) {
        cout << "This function only add empty columns at the end" << endl;
    } else if(new_size == M.dim(1)){
        return M;
    }
    VECTOR_int cc(new_size + 1, M.col_ptr(M.dim(1)));
    cc(MV_VecIndex(0,M.dim(1))) = M.col_ptr(MV_VecIndex());
    return CompCol_Mat_double(M.dim(0), new_size, M.NumNonzeros(),
            M.val(MV_VecIndex()), M.row_ind(MV_VecIndex()), cc);
}		/* -----  end of function resize_columns  ----- */


/*!
 *  \brief Concatenate a matrix and additional blocks with their number of columns
 *
 *  Take the original matrix A and add all blocks from B at the end, filling st_cols columns
 *
 *  \param A: original matrix
 *  \param B: vector of blocks to be added
 *  \param st_cols: size of the added blocks in columns
 *
 *  \return the concatenated matrix
 *
 */
    CompCol_Mat_double
concat_columns ( CompCol_Mat_double &A, std::vector<CompCol_Mat_double> &B, std::vector<int> st_cols )
{
    for(int i = 1; i < st_cols.size(); ++i)
        if(st_cols[i]<st_cols[i-1]) throw -970;

    int total_columns = A.dim(1), total_nz = A.NumNonzeros();
    int current_column = A.dim(1);

    for(int i=0; i < B.size(); i++){
        if ( st_cols[i] < current_column ) throw -971;
        total_columns += (st_cols[i] - current_column) + B[i].dim(1);
        total_nz += B[i].NumNonzeros();

        current_column = total_columns;
    }

    VECTOR_int cc(total_columns + 1), cr(total_nz);
    VECTOR_double cv(total_nz);

    cc(MV_VecIndex(0, A.dim(1))) = A.col_ptr(MV_VecIndex());
    cr(MV_VecIndex(0, A.NumNonzeros() - 1)) = A.row_ind(MV_VecIndex());
    cv(MV_VecIndex(0, A.NumNonzeros() - 1)) = A.val(MV_VecIndex());

    current_column = A.dim(1);
    int current_nz = A.NumNonzeros();
    for(int i=0; i < B.size(); i++){

        // if we add the column after a blank
        if(st_cols[i] > current_column){
            for(int j = current_column; j < st_cols[i]; j++){
                cc(j) = current_nz;
            }
            current_column = st_cols[i];
        }

        for(int j = 0; j < B[i].dim(1); j++){
            cc(j + current_column) = B[i].col_ptr(j) + current_nz;
        }
        cr(MV_VecIndex(current_nz, current_nz + B[i].NumNonzeros() -1)) = B[i].row_ind(MV_VecIndex());
        cv(MV_VecIndex(current_nz, current_nz + B[i].NumNonzeros() -1)) = B[i].val(MV_VecIndex());

        current_column += B[i].dim(1);
        current_nz += B[i].NumNonzeros();

    }
    cc(total_columns) = current_nz;

    return CompCol_Mat_double(A.dim(0), total_columns, total_nz, cv, cr, cc);
}		/* -----  end of function concat_columns  ----- */


/*!
 *  \brief Compute sparse Matrix-vector product
 *
 *  Compute sparse Matrix-vector product
 *
 *  \param M: sparse matrix
 *  \param V: sparse vector
 *
 *  \return the product sparse vector
 *
 */
    MV_ColMat_double
smv ( CompRow_Mat_double &M, MV_ColMat_double &V )
{
    assert(M.dim(1) == V.dim(0));
    //std::vector<VECTOR_double> R;
    MV_ColMat_double R(M.dim(0), V.dim(1));
    for(int k = 0; k < V.dim(1); k++){
        MV_Vector_double t = M*V(k);
        //R.push_back( M*V[k] );
        R.setCol(t, k);
    }
    return R;
}		/* -----  end of function smv  ----- */


/*!
 *  \brief Compute subpart of a sparse Matrix-vector product
 *
 *  Compute subpart of a sparse Matrix-vector product
 *
 *  \param M: sparse matrix
 *  \param col_ind: column index in the matrix to use for the product
 *  \param V: sparse vector
 *
 *  \return subpart of the product sparse vector
 *
 */
    MV_ColMat_double
spsmv ( CompRow_Mat_double &M, std::vector<int> &col_ind, MV_ColMat_double &V )
{
    MV_ColMat_double V_(M.dim(1), V.dim(1));
    for(int j = 0; j < V_.dim(1); j++){
        for(int i = 0; i < V_.dim(0); i++){
            V_(i,j) = V(col_ind[i], j);
        }
    }

    MV_ColMat_double R(M.dim(0), V.dim(1));
    for(int k = 0; k < V.dim(1); k++){
        MV_Vector_double t = M*V_(k);
        R.setCol(t, k);
    }
    return R;
}		/* -----  end of function spsmv  ----- */

/*!
 *  \brief Compute a dense Matrix-matrix product
 *
 *  Compute a dense Matrix-matrix product
 *
 *  \param L: first matrix
 *  \param R: second matrix
 *
 *  \return the dense product matrix
 *
 */
MV_ColMat_double gemmColMat(MV_ColMat_double &L, MV_ColMat_double &R)
{
    return gemmColMat(L, R, false, false);
}		/* -----  end of function gemmColMat  ----- */

/*!
 *  \brief Compute a dense Matrix-matrix product
 *
 *  Compute a dense Matrix-matrix product, where each matrix can be
 *  transposed
 *
 *  \param L: first matrix
 *  \param R: second matrix
 *  \param transL: L transposed or not
 *  \param transR: R transposed or not
 *
 *  \return the dense product matrix
 *
 */
MV_ColMat_double gemmColMat(MV_ColMat_double &L, MV_ColMat_double &R, bool transL, bool transR)
{
    assert(L.dim(1) ==  R.dim(0));

    MV_ColMat_double C(L.dim(0), R.dim(1));

    int ierr = 0;
    char no = 'N';
    //char trans = 'T';
    char tL = transL ? 'T' : 'N';
    char tR = transR ? 'T' : 'N';
    double alpha, beta;

    alpha = 1;
    beta  = 0;

    double *l_ptr = L.ptr();
    double *r_ptr = R.ptr();
    double *c_ptr = C.ptr();

    int rA = L.dim(0);
    int cA = L.dim(1);
    int rB = R.dim(0);
    int cB = R.dim(1);

    dgemm_(&tL, &tR, &rA, &cB, &cA, &alpha, l_ptr, &rA, r_ptr, &rB, &beta, c_ptr, &rA);
    return C;
}		/* -----  end of function gemmColMat  ----- */

/*!
 *  \brief Returns the upper triangular matrix
 *
 *  Returns the upper triangular matrix
 *
 *  \param M: the complete matrix
 *
 *  \return the upper triangular matrix
 *
 */
MV_ColMat_double upperMat(MV_ColMat_double &M){
    MV_ColMat_double R = M;
    for(int i = 0; i < M.dim(0); i++)
        for(int j = 0; j < i; j++)
            R(i,j) = 0;
    return R;
}		/* -----  end of function upperMat  ----- */


/*!
 *  \brief Sparse Matrix-Matrix product
 *
 *  Returns the sparse Matrix-Matrix product of 2 sparse matrices
 *
 *  \param A: the first matrix
 *  \param BT: the second matrix
 *
 *  \return the sparse Matrix-Matrix product
 *
 */
CompRow_Mat_double spmm_overlap ( CompRow_Mat_double &A, CompRow_Mat_double &BT ) {
    MKL_INT p, jp, j, kp, k, i, nz = 0, anz, *Cp, *Cj, *Ap, *Aj, *Bp, m, n,
    bnz, *Bj;

    double t1 = MPI_Wtime();
    double *Ax, *Bx, *Cx;

    m = A.dim(0);
    Ap = A.rowptr_ptr();
    Aj =  A.colind_ptr();
    Ax =  A.val_ptr();
    anz = A.NumNonzeros();

    n = BT.dim(0);
    Bp = BT.rowptr_ptr();;
    Bj = BT.colind_ptr();
    Bx = BT.val_ptr();
    bnz = BT.NumNonzeros();

    MKL_INT *xb = (int*)calloc(m,sizeof(MKL_INT));
    double *x = (double*)calloc(m,sizeof(double));
    MKL_INT nzmax = (anz + bnz) * 2;
    MKL_INT nzt = (anz + bnz);

    Cp = (MKL_INT*) calloc(n+1,sizeof(MKL_INT));
    Cj = (MKL_INT*) calloc(nzmax,sizeof(MKL_INT));
    Cx = (double*) calloc(nzmax,sizeof(double));


    //~ #pragma omp parallel for firstprivate(xb,x,k,j,p) private(i)
    for (i = 0; i < m; i++) {
        if ( ( (nz + n) > nzmax ) ) {
            nzmax += nzt ;
            cout << "nzmax " << nzmax << " "<< nzt <<endl;
            Cj = (MKL_INT*)realloc(Cj, sizeof(MKL_INT)*nzmax);
            Cx = (double*)realloc(Cx, sizeof(double)*nzmax);
            if(Cj == NULL || Cx == NULL){
                cout << "out of Memory " << nzmax << " "<< anz <<endl;
                exit(0);
            }
        }

        Cp[i] = nz; /* row i of C starts here */
        for (int jp = Ap[i]; jp < Ap[i + 1]; jp++) {
            j = Aj[jp];
            for (int kp = Bp[j]; kp < Bp[j + 1]; kp++) {
                k = Bj[kp]; /* B(i,j) is nonzero */
                if(k != i) {
                    if (xb[k] != i +1) {
                        xb[k] = i +1; /* i is new entry in column j */
                        Cj[nz++] = k; /* add i to pattern of C(:,j) */
                        x[k] = Ax[jp] * Bx[kp]; /* x(i) = beta*A(i,j) */
                    } else {
                        x[k] += (Ax[jp] * Bx[kp]); /* i exists in C(:,j) already */
                    }
                }
            }
        }
        for (p = Cp[i]; p < nz; p++){
            Cx[p] = abs(x[Cj[p]]);
        }
    }

    cout << "SPMM Multiplication done. Time spent: " <<  t1 - MPI_Wtime() << endl;
    Cp[m] = nz; /* finalize the last row of C */
    Cj = (MKL_INT*)realloc(Cj, sizeof(MKL_INT)*nz);
    Cx = (double*)realloc(Cx, sizeof(double)*nz);

    CompRow_Mat_double C (n , n, nz , Cx, Cp, Cj);

    delete Cp, Cj,Cx;
    delete xb, x;

    return C;
}       /* -----  end of function spmm  ----- */
