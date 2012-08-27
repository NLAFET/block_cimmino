#include "abcd.h"

///@TODO Add the case of manual partitioning
int abcd::partition()
{
    unsigned nbrows_per_part;
    unsigned row_sum = 0;

    if (nbparts == 0)
        return -13;

    nbrows_per_part = (unsigned)(m / nbparts);

    if (partitioning_type == 2) {
        strow = ArrayXi(nbparts);
        nbrows = ArrayXi(nbparts);

        nbrows.setConstant(nbrows_per_part);

        nbrows(nbparts - 1) += m - nbrows.sum();
        for (unsigned k = 0; k < nbparts; k++){
            strow(k) = row_sum;
            row_sum += nbrows(k);
        }
    }

    return 0;
}
