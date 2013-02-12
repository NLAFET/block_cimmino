/*
 * =====================================================================================
 *
 *       Filename:  solve.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/06/2013 11:36:33 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */


#include <abcd.h>


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::solveABCD
 *  Description:  
 * =====================================================================================
 */
    void
abcd::solveABCD ( MV_ColMat_double &b )
{
    // if not created yet, do it!
    if( S.dim(0) == 0 ){
        S = abcd::buildS();
    }
}		/* -----  end of function abcd::solveABCD  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::buildS
 *  Description:  
 * =====================================================================================
 */
    Coord_Mat_double
abcd::buildS (  )
{
    Coord_Mat_double shur;
    std::vector<int> sr, sc;
    std::vector<double> sv;

    Xk = MV_ColMat_double(n, 1, 0);
    MV_ColMat_double b(m, 1, 0); 

    glob_to_local[n_o + 1]; // just in case it does not exist.

    for( int i = 0; i < size_c; i++){

        // set xk = [0 I] 
        if(glob_to_local[n_o + i] != 0){
            Xk(glob_to_local[n_o + i - 1], 0) = 0;
            Xk(glob_to_local[n_o + i], 0) = 1;
        }

        // solve with b = 0; and xk_0 = xk
        bcg(b);

        if(inter_comm.rank() == 0){
            VECTOR_double xx(size_c, 0);
            for(int k = 0; k < partitions.size(); k++) {
                for( int j = 0; j < size_c; j++){
                    if(glob_to_local[n_o+j]!=0){
                        //cout << n_o + j << " "  << glob_to_local[n_o + j] << " " <<Xk(glob_to_local[n_o + j], 0) << endl;
                    }
                }
                exit(0);
            }
        }
    }

    return shur;
}		/* -----  end of function abcd::buildS  ----- */
