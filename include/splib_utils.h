/*!
 * \file splib_utils.h
 * \brief Header for utils functions added to Sparselib
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

double squaredNorm(CompRow_Mat_double &M);
double squaredNorm(VECTOR_double &V);
double squaredNorm(VECTOR_double &V, VECTOR_int &I);
double infNorm(VECTOR_double &V);
double infNorm(MV_ColMat_double &V);
double infNorm(Coord_Mat_double &M);
CompRow_Mat_double CSC_extractByIndices  (CompRow_Mat_double &M, std::vector<int> row_index);
CompCol_Mat_double CSC_middleRows (CompRow_Mat_double &M, int st_row, int nb_rows);
CompCol_Mat_double sub_matrix (CompCol_Mat_double &M, std::vector<int> &ci);
VECTOR_double middleCol(CompCol_Mat_double &M, int col_num);
VECTOR_double middleCol(CompCol_Mat_double &M, int col_num, VECTOR_int &ind);
CompRow_Mat_double smmtm (CompRow_Mat_double &A, CompRow_Mat_double &B, double prune);
CompRow_Mat_double smmtm (CompCol_Mat_double &A, CompCol_Mat_double &B, double prune);
CompRow_Mat_double smmtm (CompRow_Mat_double &A, CompRow_Mat_double &B);
CompRow_Mat_double smmtm (CompCol_Mat_double &A, CompCol_Mat_double &B);
CompRow_Mat_double spmm (CompRow_Mat_double &A, CompRow_Mat_double &B);
CompCol_Mat_double csc_transpose ( CompCol_Mat_double &M );
CompCol_Mat_double csc_transpose ( CompRow_Mat_double &M );
CompRow_Mat_double csr_transpose ( CompCol_Mat_double &M );
CompRow_Mat_double csr_transpose ( CompRow_Mat_double &M );
CompCol_Mat_double resize_columns ( CompCol_Mat_double &M, int new_size );
CompCol_Mat_double concat_columns ( CompCol_Mat_double &A, std::vector<CompCol_Mat_double> &B, std::vector<int> st_cols );
MV_ColMat_double smv ( CompRow_Mat_double &M, MV_ColMat_double &V );
MV_ColMat_double spsmv ( CompRow_Mat_double &M, std::vector<int> &col_ind, MV_ColMat_double &V );
MV_ColMat_double gemmColMat(MV_ColMat_double &L, MV_ColMat_double &R);
MV_ColMat_double gemmColMat(MV_ColMat_double &L, MV_ColMat_double &R, bool transL, bool transR);
MV_ColMat_double upperMat(MV_ColMat_double &M);
CompRow_Mat_double spmm_overlap ( CompRow_Mat_double &A, CompRow_Mat_double &BT );
