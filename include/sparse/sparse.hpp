#ifndef _SPARSE_HPP_
#define _SPARSE_HPP_

#include <vector>
#include <iostream>

//! Sparse matrix class
/*!
 * @tparam ValueType the type of the entries of the matrix, for the moment only double is supported
 * @tparam IndicesType the type of vector that stores row/column indices, default is size_t
 * @tparam base the base of the vector indices, 0 for c-based, 1 for fortran-based
 */
template< typename ValueType, typename IndicesType, int base = 0>
class SparseMat {
    protected:
        ValueType* val_; /*!< data values (nz_ elements)*/
        IndicesType* rows_; /*!< row indices or pointers, depends on the storage format (m_ or nz_ elements) */
        IndicesType* cols_; /*!< column indices or pointers, depends on the storage format (n_ or nz_ elements) */

        IndicesType nz_; /*!< the number of nonzero elements */
        IndicesType m_; /*!< the number of rows */
        IndicesType n_; /*!< the number of columns */

        bool isSymmetric; /*!< defines wether the matrix is symmetric or not */
        bool isReference; /*!< defines if the matrix is a reference to dynamic arrays */

    public:
        /*!
         * Basic constructor, initialize everything to zero
         */
        SparseMat() : isSymmetric(false), m_(0), n_(0), nz_(0) { }

        //! Returns a pointer to the first element of the values vector
        ValueType *pVal() const { return &val_[0]; }
        //! Returns a pointer to the first element of the rows vector
        IndicesType *pRows() const { return &rows_[0]; }
        //! Returns a pointer to the first element of the cols vector
        IndicesType *pCols() const { return &cols_[0]; }

        //! Return the number of rows
        IndicesType nRows() const { return m_; };
        //! Return the number of columns
        IndicesType nCols() const { return n_; };
        //! Return the number of nonzero elements
        IndicesType nnz() const { return nz_; };

        //! Get (read-only) the element at the i-th row and j-th column
        /*!
         * @param i the row indice
         * @param j the column indice
         * @return the element at the i-th row and j-th column
         */
        ValueType get(IndicesType i, IndicesType j) const;

        //! Access (read/write) the element at the i-th row and j-th column
        /*!
         * @param i the row indice
         * @param j the column indice
         * @return the element at the i-th row and j-th column
         */
        ValueType& operator() (IndicesType i, IndicesType j) const;
};

#endif // _SPARSE_HPP_
