/*
 * =====================================================================================
 *
 *       Filename:  splib_utils.cpp
 *
 *    Description:  Some additions to SparseLib++
 *
 *        Version:  1.0
 *        Created:  01/23/2013 06:47:59 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include <abcd.h>
#include <Eigen/src/misc/blas.h>

//extern "C"
//{
    //void amub_ (int *nrow,int *ncol,int *job,double *a,int *ja,int *ia,double *b,int *jb,int *ib,
           //double *c,int *jc,int *ic,int *nzmax,int *iw,int *ierr );
//}

/// Sum nonzeros from parts
int sum_nnz(int res, CompRow_Mat_double M)
{
    return res += M.NumNonzeros();
}
int sum_rows(int res, CompRow_Mat_double M)
{
    return res += M.dim(0);
}
int sum_cols(int res, CompRow_Mat_double M)
{
    return res += M.dim(0);
}
//bool comp_cols(Eigen::SparseMatrix<double, RowMajor> L, Eigen::SparseMatrix<double, RowMajor> R)
//{
    //return (L.cols() < R.cols());
//}
//
double infNorm(Coord_Mat_double &M){
    double max = 0;
    double *pt = M.val_ptr();
    double v;
    for (int i = 0; i < M.NumNonzeros(); i++){
        //if(abs(V(i)) >= max) max = abs(V(i));
        v = abs(pt[i]);
        max = v > max ? v : max;
    }
    return max;
}
double infNorm(VECTOR_double &V){
    double max = 0;
    for (int i = 0; i < V.size(); i++){
        //if(abs(V(i)) >= max) max = abs(V(i));
        max = abs(V(i)) > max ? abs(V(i)) : max;
    }
    return max;
}

double infNorm(MV_ColMat_double &V){
    double max = 0;
    double *v_ptr = V.ptr();
    for (int i = 0; i < V.dim(0)*V.dim(1); i++){
        if(abs(v_ptr[i]) >= max) max = abs(v_ptr[i]);
    }
    return max;
}

double squaredNorm(VECTOR_double &V){
    double sum = 0;
    double e;
    for (int i = 0; i < V.size(); i++){
        e = V(i);
        sum += e*e;
    }
    return sum;
}

double squaredNorm(VECTOR_double &V, VECTOR_int &I){
    double sum = 0;
    double e;
    for (int i = 0; i < I.size(); i++){
        e = V(I(i));
        sum += e*e;
    }
    return sum;
}

double squaredNorm(CompRow_Mat_double &M)
{
    double sum = 0;
    for (int i = 0; i < M.NumNonzeros(); i++){
        double v = M.val(i);
        sum += v*v;
    }
    return sum;
}

CompRow_Mat_double CSR_middleRows (CompRow_Mat_double &M, int st_row, int nb_rows) {
    int st_index, ed_index;

    // starting index in JCN
    st_index = M.row_ptr(st_row);
    // last index in JCN
    ed_index = M.row_ptr(st_row + nb_rows) - 1;

    VECTOR_int sub_row_vect(nb_rows + 1);
    sub_row_vect() = M.row_ptr(MV_VecIndex(st_row, st_row + nb_rows));

    int starting_point = M.row_ptr(st_row);
    VECTOR_double vv = M.val(MV_VecIndex(st_index, ed_index));


    for(int i = 0; i <= nb_rows; i++){
        sub_row_vect[i] -= starting_point;
    }

    VECTOR_int sub_col_vect(ed_index - st_index + 1);
    sub_col_vect() = M.col_ind(MV_VecIndex(st_index, ed_index));

    VECTOR_double sub_val_vect(ed_index - st_index + 1);
    sub_val_vect() = M.val(MV_VecIndex(st_index, ed_index));

    return CompRow_Mat_double( nb_rows, M.dim(1), ed_index - st_index + 1,
            sub_val_vect,
            sub_row_vect,
            sub_col_vect
            );
}

CompCol_Mat_double CSC_middleRows (CompRow_Mat_double &M, int st_row, int nb_rows) {
    return CompCol_Mat_double(CSR_middleRows(M, st_row, nb_rows));
}

VECTOR_double middleCol(CompCol_Mat_double &M, int col_num){
    VECTOR_int ind;
    return middleCol(M, col_num, ind);

}

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
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sub_matrix
 *  Description:  
 * =====================================================================================
 */
    CompCol_Mat_double
sub_matrix ( CompCol_Mat_double &M, std::vector<int> &ci )
{
    CompCol_Mat_double SM;
    int st_col, ed_col;
    int c = 0, j = 0, nzc;

    VECTOR_int sm_r;//(M.dim(0)*ci.size());
    VECTOR_int sm_c(ci.size() + 1, 0);
    VECTOR_double sm_v;//(M.dim(0)*ci.size());

    std::vector<int> v_sm_r;
    std::vector<double> v_sm_v;

    for(int k=0; k < ci.size(); k++){
        st_col = M.col_ptr(ci[k]);
        ed_col = M.col_ptr(ci[k] + 1);

        nzc = ed_col - st_col;

        // if the column is empty
        if (st_col == ed_col) continue;

        c += nzc;
        j++;
    }

    v_sm_r.reserve(nzc);
    v_sm_v.reserve(nzc);

    c = 0; j = 0;


    // for all columns in ci
    for(int k=0; k < ci.size(); k++){
        st_col = M.col_ptr(ci[k]);
        ed_col = M.col_ptr(ci[k] + 1);

        nzc = ed_col - st_col;

        // if the column is empty
        if (st_col == ed_col) continue;

        sm_c[j] = c;
        //sm_r(MV_VecIndex(c, c + nzc - 1)) = M.row_ind(MV_VecIndex(st_col, ed_col - 1));
        //sm_v(MV_VecIndex(c, c + nzc - 1)) = M.val(MV_VecIndex(st_col, ed_col - 1));
        //
        int pos = st_col;
        for(int i = c; i < c + nzc; i++){
            v_sm_r.push_back(M.row_ind(pos));
            v_sm_v.push_back(M.val(pos));
            pos++;
        }

        c += nzc;
        j++;
    }
    if (c==0) return SM;
    sm_c[j] = c;
    SM = CompCol_Mat_double(M.dim(0), j, c,
            //sm_v(MV_VecIndex(0,c-1)), sm_r(MV_VecIndex(0,c-1)), sm_c);
            &v_sm_v[0], &v_sm_r[0], sm_c.ptr());

    return SM;
}		/* -----  end of function subMatrix  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  smmtm
 *  Description:  Sparse Matrix-Matrix^T Multiplication
 *                C = A * B^T
 * =====================================================================================
 */
CompRow_Mat_double smmtm (CompRow_Mat_double &A, CompRow_Mat_double &B)
{
    CompRow_Mat_double BT = csr_transpose(B);

    int nzmax = A.dim(0) * B.dim(0);

    VECTOR_int iw(BT.dim(1), 0);
    std::map<int,int> jc, ic;
    std::map<int,double> c;
    ic[0] = 0;

    double scal;
    int i, ka, kb, j, k, jcol, jpos, len=0;

    if(A.dim(1) != BT.dim(0)) throw -980;

    for(i = 0; i < A.dim(0); i++){
        // ith row
        for(ka = A.row_ptr(i); ka < A.row_ptr(i+1); ka++){

            scal = A.val(ka);
            j = A.col_ind(ka);

            for(kb = BT.row_ptr(j); kb < BT.row_ptr(j+1); kb++){
                jcol = BT.col_ind(kb);
                jpos = iw(jcol);

                if(jpos == 0){
                    len = len + 1;
                    if (nzmax < len) {
                        throw -981;
                    }

                    jc[len - 1] = jcol;
                    iw[jcol] = len;

                    c[len - 1] = scal*BT.val(kb);
                    //cout << "->" << len << "\t" << c[len -1 ] << "\t" << jc[len -1] <<  endl;
                } else {
                    c[jpos - 1] = c[jpos - 1] + scal * BT.val(kb);
                    //cout << "--" << jpos << "\t" << c[jpos -1 ] << "\t" << jc[jpos -1] <<  endl;

                    //cout << i+1 << " " << jcol+1 << "\t" << c[jpos -1 ] << endl;
                }
            }

        }
        int l = len;
        for(k = ic[i]; k < l; k++){
            if(c[k] == 0){
                jc.erase(k);
                c.erase(k);
                len--;
            }
            iw(jc[k]) = 0;
        }

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

}

CompRow_Mat_double smmtm (CompCol_Mat_double &A, CompCol_Mat_double &B)
{
    CompRow_Mat_double Ar = A;
    CompRow_Mat_double Br = B;
    return smmtm(Ar, Br);
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  spmm
 *  Description:  
 * =====================================================================================
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
    int i, ka, kb, j, k, jcol, jpos, len=0;

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
}		/* -----  end of function spmm  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  csc_transpose
 *  Description:  
 * =====================================================================================
 */
    CompCol_Mat_double
csc_transpose ( CompCol_Mat_double &M )
{
    CompRow_Mat_double Mt = csr_transpose(M);
    return CompCol_Mat_double(Mt);

}		/* -----  end of function csc_transpose  ----- */
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

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  csr_transpose
 *  Description:  
 * =====================================================================================
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

    CompRow_Mat_double
csr_transpose ( CompRow_Mat_double &M )
{
    return CompRow_Mat_double(csc_transpose(M));
}		/* -----  end of function csr_transpose  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  resize_columns
 *  Description:  
 * =====================================================================================
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


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  concat_columns
 *  Description:  
 * =====================================================================================
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


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  smv
 *  Description:  
 * =====================================================================================
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
}		/* -----  end of function smdm  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  spsmv
 *  Description:  
 * =====================================================================================
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

MV_ColMat_double gemmColMat(MV_ColMat_double &L, MV_ColMat_double &R)
{
    return gemmColMat(L, R, false, false);
}

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
}

MV_ColMat_double upperMat(MV_ColMat_double &M){
    MV_ColMat_double R = M;
    for(int i = 0; i < M.dim(0); i++)
        for(int j = 0; j < i; j++)
            R(i,j) = 0;
    return R;
}
