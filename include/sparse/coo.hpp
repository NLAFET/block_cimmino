#ifndef _COO_HPP_
#define _COO_HPP_

#include <vector>
#include <iostream>
#include <stdlib.h>

#include "sparse.hpp"

using namespace std;

//! Coordinate storage format sparse matrix
/*!
 * @tparam ValueType the type of the entries of the matrix, for the moment only double is supported
 * @tparam IndicesType the type of vector that stores row/column indices, default is size_t
 * @tparam base the base of the vector indices, 0 for c-based, 1 for fortran-based
 */
template<
    typename ValueType = double,
    typename IndicesType = int,
    int base = 0
    >
class CooMatrix : public SparseMat<ValueType,IndicesType,base> {
    public:
        CooMatrix( IndicesType m,  IndicesType n,  IndicesType nz,
                 IndicesType* irn,  IndicesType* jcn,  ValueType* val,
                 bool reference)
        {
            this->m_ = m;
            this->n_ = n;
            this->nz_ = nz;
            if (reference) {
                this->val_  = val;
                this->rows_ = irn;
                this->cols_ = jcn;
            } else {
                this->val_ = new ValueType[this->nz_];
                this->rows_ = new IndicesType[this->nz_];
                this->cols_ = new IndicesType[this->nz_];
            }
        }
        CooMatrix( IndicesType m,  IndicesType n,  IndicesType nz,
                 IndicesType* irn,  IndicesType* jcn,  ValueType* val)
        {
            CooMatrix(m, n, nz, irn, jcn, val, false);
        }


        //! Returns the entry at position @i
        ValueType val(IndicesType i) const { return this->val_[i]; }
        //! Returns the row indice at position @i
        IndicesType row(IndicesType i) const { return this->rows_[i]; }
        //! Returns the column indice at position @i
        IndicesType col(IndicesType i) const { return this->cols_[i]; }

        //! Slow get (read-only) of the element at the i-th row and j-th column
        /*!
         * @param i the row indice
         * @param j the column indice
         * @return the element at the i-th row and j-th column
         */
        ValueType get(IndicesType i, IndicesType j) const;

        //! Slow access (read/write) to the element at the i-th row and j-th column
        /*!
         * @param i the row indice
         * @param j the column indice
         * @return the element at the i-th row and j-th column
         */
        ValueType& operator() (IndicesType i, IndicesType j);

        //! Print all of the matrix entries and uses base = 1
        void show();

        //! Print all of the matrix entries and uses @customBase as base
        void show(int customBase);
};

template<typename ValueType, typename IndicesType, int base>
ValueType& CooMatrix<ValueType,IndicesType,base>::operator() (IndicesType i, IndicesType j)
{
    for (int t = 0; t < this->nz_; t++)
        if (this->rows_[t] == i && this->cols_[t] == j) return this->val_[t];

    std::cerr << "Array element (" << i << "," << j ;
    std::cerr << ") not in sparse structure -- cannot access." << "\n";
    /// @todo: define exceptions
    exit(1);
    return this->val_[0];
}

template<typename ValueType, typename IndicesType, int base>
ValueType CooMatrix<ValueType,IndicesType,base>::get(IndicesType i, IndicesType j) const
{
    for (int t = 0; t < this->nz_; t++)
        if (this->rows_[t] == i && this->cols_[t] == j) return this->val_[t];

    // default behavior if we do not have this entry
    if (i < this->m_ && j < this->n_) return 0.0;

    return this->val_[0];
}

template<typename ValueType, typename IndicesType, int base>
ostream& operator << (ostream & os, const CooMatrix<ValueType,IndicesType,base> & mat)
{
         IndicesType nnz = mat.nnz();
         IndicesType M = mat.nRows();
         IndicesType N = mat.nCols();
         IndicesType rowp1, colp1;
         int flag = 0;


        std::ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
        std::ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);

         int oldp = os.precision(12);

//       Loop through Nonzeros
         for (IndicesType j = 0; j < nnz ; j++) 
         {
            rowp1 = mat.row(j) +1;
            colp1 = mat.col(j) +1;
            if ( rowp1 == M && colp1 == N ) flag = 1;
            os.width(14);
            os <<  rowp1 ; os << "    " ;
            os.width(14);
            os <<  colp1 ; os << "    " ;
            os.width(20);
            os <<  mat.val(j) << "\n";
         }
         if (flag == 0) 
         {
            os.width(14);
            os <<  M ; os << "    " ;
            os.width(14);
            os <<  N ; os << "    " ;
            os.width(20);
            os <<  mat.get(M-1,N-1) << "\n";
         }
         os.setf(olda,ios::adjustfield);
         os.setf(oldf,ios::floatfield);
         os.precision(oldp);
         return os;
}

#endif // _COO_HPP_
