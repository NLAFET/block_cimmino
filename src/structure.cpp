#include "abcd.h"

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

    if(partitioning_type == 2) {
        strow = ArrayXi(nbparts);
        nbrows = ArrayXi(nbparts);

        nbrows.setConstant(nbrows_per_part);

        nbrows(nbparts - 1) += m - nbrows.sum();
        for(unsigned k = 0; k < nbparts; k++) {
            strow(k) = row_sum;
            row_sum += nbrows(k);
        }
    }

}

void abcd::analyseFrame()
{
    for(unsigned k = 0; k < nbparts; k++) {
        // Our k-th partition
        SparseMatrix<double, ColMajor> part = mtx.middleRows(strow[k], nbrows[k]);

        double t1, t2;
        // Build the column index of part
        std::vector<int> ci;
        int j = 0;
        for(int i = 1; i <= part.outerSize(); i++) {
            if(part.outerIndexPtr()[i] != part.outerIndexPtr()[i - 1])
                ci.push_back(j);
            j++;
        }
        column_index.push_back(ci);

        int *last = std::unique(part.outerIndexPtr(), part.outerIndexPtr() + part.outerSize() + 1);
        parts.push_back(SparseMatrix<double, RowMajor>(part.middleCols(0, ci.size())));
    }
}