#include "abcd.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::augmentMatrix
 *  Description:  Augments the matrix and build the C part in [A C]
 * =====================================================================================
 */
    void
abcd::augmentMatrix ( std::vector<CompCol_Mat_double> &M)
{
    /*
     * Which augmentation to use:
     */
    if(icntl[Controls::aug_type] == 0){
        //[> No augmentation <]
        return;
    } else if (icntl[Controls::aug_type] == 1){
        /*
         * C_ij/-I augmentation
         */
        cijAugmentMatrix(M);
    } else if (icntl[Controls::aug_type] == 2){
        /*
         * A_ij/-A_ji augmentation
         */
        aijAugmentMatrix(M);

    } else {
        throw std::runtime_error("Unkown augmentation scheme.");
    }

}// [> -----  end of function abcd::augmentMatrix  ----- <]
