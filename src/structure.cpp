#include <abcd.h>

///@TODO Add the case of manual partitioning
void abcd::partitionMatrix()
{
    unsigned nbrows_per_part;
    unsigned row_sum = 0;

    if(nbparts == 0)
        throw - 13;
    if(nbparts < parallel_cg)
        throw - 101;

    nbrows_per_part = ceil(float(m)/float(nbparts));

    switch(partitioning_type){
        
        /*-----------------------------------------------------------------------------
         *  Uniform partitioning with a given nbrows
         *-----------------------------------------------------------------------------*/
        case 1:
            strow = ArrayXi(nbparts);

            for(unsigned k = 0; k < nbparts; k++) {
                strow(k) = row_sum;
                row_sum += nbrows(k);
            }
            break;

        /*-----------------------------------------------------------------------------
         *  Uniform partitioning with only nbparts as input (generates nbrows)
         *-----------------------------------------------------------------------------*/
        case 2:
            strow = ArrayXi(nbparts);
            nbrows = ArrayXi(nbparts);

            nbrows.setConstant(nbrows_per_part);

            nbrows(nbparts - 1) += m - nbrows.sum();
            for(unsigned k = 0; k < nbparts; k++) {
                strow(k) = row_sum;
                row_sum += nbrows(k);
            }
            break;
        /*-----------------------------------------------------------------------------
         *  PaToH partitioning
         *-----------------------------------------------------------------------------*/
        //case 3:

    }

}

void abcd::analyseFrame()
{

    std::vector<Eigen::SparseMatrix<double, ColMajor> > loc_parts;
    std::vector<int> ci_sizes;
    for(unsigned k = 0; k < nbparts; k++) {
        // Our k-th partition
        SparseMatrix<double, ColMajor> part = mtx.middleRows(strow[k], nbrows[k]);
        loc_parts.push_back(part);
    }


    if(use_abcd){
        abcd::augmentMatrix(loc_parts);
    }

    for(unsigned k = 0; k < nbparts; k++) {
        double t1, t2;
        // Build the column index of part
        std::vector<int> ci;
        int j = 0;
        for(int i = 1; i <= loc_parts[k].outerSize(); i++) {
            if(loc_parts[k].outerIndexPtr()[i] != loc_parts[k].outerIndexPtr()[i - 1])
                ci.push_back(j);
            j++;
        }
        column_index.push_back(ci);
        ci_sizes.push_back(ci.size());

        //int *last = std::unique(part.outerIndexPtr(), part.outerIndexPtr() + part.outerSize() + 1);
        //parts.push_back(SparseMatrix<double, RowMajor>(part.middleCols(0, ci.size())));
        int *last = std::unique(loc_parts[k].outerIndexPtr(),
                loc_parts[k].outerIndexPtr() + loc_parts[k].outerSize() + 1);
        parts.push_back(SparseMatrix<double, RowMajor>(loc_parts[k].middleCols(0, ci.size())));
        //loc_parts[k] = NULL;
    }

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::augmentMatrix
 *  Description:  Augments the matrix and build the C part in [A C]
 * =====================================================================================
 */
    void
abcd::augmentMatrix ( std::vector<Eigen::SparseMatrix<double, ColMajor> > &loc_parts)
{
    /*
     * Which augmentation to use:
     */
    if(icntl[10] == 0){
        /* No augmentation */
        return;
    } else if (icntl[10] == 1){
        /*
         * C_ij/-I augmentation
         */
        
        // build C_ij

    } else if (icntl[10] == 2){
        /*
         * A_ij/-A_ji augmentation
         */

    } else if (icntl[10] == 3){
        /*
         * SVD augmentation
         */

    }

}		/* -----  end of function abcd::augmentMatrix  ----- */
