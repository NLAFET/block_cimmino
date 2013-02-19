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
    if(inter_comm.rank() == 0)
        cout << " [-] Building S" << endl;


    glob_to_local[n_o + 1]; // just in case it does not exist.

    //std::vector<int> indices;
    //boost::copy(glob_to_local | map_keys, std::back_inserter(indices));

    bool use_xk_t = use_xk;

    use_xk = true;
    //
    for( int i = 0; i < 1; i++){

        // set xk = [0 I] 
        Xk = MV_ColMat_double(n, 1, 0);
        // Retrieve all keys
        //cout << inter_comm.rank() << " " << indices[0] << " " << indices.back() << " " << n_o + i << endl;
        std::map<int,int>::iterator it = glob_to_local.find(n_o + i);

        if(it!=glob_to_local.end()){
            Xk(glob_to_local[n_o + i], 0) = 1;
        }

        MV_ColMat_double b(m, 1, 0); 

        //if(it!=glob_to_local.end()){
            //int st = 0;
            //for(int p = 0; p < partitions.size(); p++){
                    //b(MV_VecIndex(st, st + partitions[p].dim(0) - 1),
                            //MV_VecIndex(0, 0)) = spsmv(partitions[p], column_index[p], Xk);
                //st += partitions[p].dim(0);
            //}
        //}

        // solve with b = 0; and xk_0 = xk
        if(inter_comm.rank() == 0)
            cout << "     Column " << i << endl;

        if(dcntl[10] == 0){

        } else {
            bcg(b);
            //if(inter_comm.rank() == 0)
                //cout << Xk << endl;
        }
        //exit(0);

        //if(inter_comm.rank() == 0){
            //VECTOR_double xx(size_c, 0);
            //for(int k = 0; k < partitions.size(); k++) {
                //for( int j = 0; j < size_c; j++){
                    //if(glob_to_local[n_o+j]!=0){
                        ////cout << n_o + j << " "  << glob_to_local[n_o + j] << " " <<Xk(glob_to_local[n_o + j], 0) << endl;
                    //}
                //}
                ////exit(0);
            //}
        //}
    }
    use_xk = false;

    return shur;
}		/* -----  end of function abcd::buildS  ----- */
