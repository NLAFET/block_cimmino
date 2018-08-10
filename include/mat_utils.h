#ifndef MAT_UTILS_HXX_
#define MAT_UTILS_HXX_

/*!
 * \file mat_utils.h
 * \brief Implementation of utils on matrix computation
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include<vector>

/*!
 *  \brief Transform a CSR matrix to CSC format
 *
 *  Transform a CSR matrix to CSC format
 *
 *  \param n_row: number of rows
 *  \param n_col: number of columns
 *  \param Ap: array of rows
 *  \param Aj: array of columns
 *  \param Ax: array of values
 *  \param Bp: array of rows of output
 *  \param Bj: array of columns of output
 *  \param Bx: array of values of output
 *
 */
template <class I, class T>
void csr_tocsc(const I n_row,
               const I n_col,
               const I Ap[],
               const I Aj[],
               const T Ax[],
               I Bp[],
               I Bi[],
               T Bx[])
{
    const I nnz = Ap[n_row];

    //compute number of non-zero entries per column of A
    std::fill(Bp, Bp + n_col, 0);

    for (I n = 0; n < nnz; n++){
        Bp[Aj[n]]++;
    }

    //cumsum the nnz per column to get Bp[]
    for(I col = 0, cumsum = 0; col < n_col; col++){
        I temp = Bp[col];
        Bp[col] = cumsum;
        cumsum += temp;
    }
    Bp[n_col] = nnz;

    for(I row = 0; row < n_row; row++){
        for(I jj = Ap[row]; jj < Ap[row+1]; jj++){
            I col = Aj[jj];
            I dest = Bp[col];

            Bi[dest] = row;
            Bx[dest] = Ax[jj];

            Bp[col]++;
        }
    }

    for(I col = 0, last = 0; col <= n_col; col++){
        I temp = Bp[col];
        Bp[col] = last;
        last = temp;
    }
} 
#endif // MAT_UTILS_HXX_

