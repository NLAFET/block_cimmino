/*
 * =====================================================================================
 *
 *       Filename:  s_utils.cpp
 *
 *    Description:  Contain the different methods related to building S and
 *    computing columns of S
 *
 *        Version:  1.0
 *        Created:  03/11/2013 02:11:05 PM
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
 *         Name:  abcd::buildS
 *  Description:  Builds a sparse S
 * =====================================================================================
 */
    Coord_Mat_double
abcd::buildS (  )
{
    Coord_Mat_double shur;
    //MV_ColMat_double S(size_c, size_c, 0);
    std::vector<int> sr, sc;
    std::vector<double> sv;
    if(inter_comm.rank() == 0)
        cout << " [->] Building S" << endl;


    std::vector<int> vc, vr;
    std::vector<double> vv;

    std::vector<int> my_cols;


    if(dcntl[10] == 0){
        for( int i = 0; i < size_c; i++){
            std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);
            if(iti!=glob_to_local.end()) my_cols.push_back(i);
        }

        std::vector<int>::iterator pos = my_cols.begin();
        std::vector<int>::iterator end_pos;

        int share = 16;
        while(pos != my_cols.end()){
            if(pos + share < my_cols.end()) end_pos = pos + share;
            else end_pos = my_cols.end();

            std::vector<int> cur_cols;

            std::copy(pos, end_pos, std::back_inserter(cur_cols));

            int mumps_share = share > 32 ? share/2 : 16;
            setMumpsIcntl(27, mumps_share);
            MV_ColMat_double sp = spSimpleProject(cur_cols);

            for( int j = 0; j < cur_cols.size(); j++){
                int c = cur_cols[j];
                for( int i = 0; i < size_c; i++){
                    if(sp(i,j)!=0){
                        vc.push_back(c);
                        vr.push_back(i);
                        vv.push_back(sp(i, j));
                    }
                }
            }
            pos = end_pos;
        }


    } else {
        for( int i = 0; i < size_c; i++){
            std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);
            if(iti!=glob_to_local.end()) my_cols.push_back(i);
        }

        //for( int j = 0; j < my_cols.size(); j++){
            //int i = my_cols[j];
        for( int i = 0; i < size_c; i++){
        //
            int my_bro = -1;
            if(inter_comm.rank() == 0) cout << " Column " << i << " out of " << size_c << endl;

            //for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
                    //it != col_interconnections.end(); it++) {
                //int k = it->second.size() - 1;
                //while( k > -1 && my_bro == -1){
                    //if(it->second[k] == glob_to_local[n_o + i]){
                        //my_bro = it->first;
                    //}
                    //k--;
                //}
            //}


            // set xk = [0 I] 
            block_size = 1;
            Xk = MV_ColMat_double(n, 1, 0);
            MV_ColMat_double b(m, 1, 0); 

            std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);
            if(iti!=glob_to_local.end()){
                Xk(glob_to_local[n_o + i], 0) = 1;

                //int st = 0;
                //for(int p = 0; p < partitions.size(); p++){
                    //MV_ColMat_double sp = spsmv(partitions[p], local_column_index[p], Xk);
                    //int pos = 0;
                    //for(int k = st; k < st + partitions[p].dim(0); k++){
                        //b(k, 0) = sp(pos, 0);
                        //pos++;
                    //}
                    //st += partitions[p].dim(0);
                //}
            }
            //MV_ColMat_double sp = Xk - coupleSumProject(0e0, b, 1e0, Xk, my_bro);

            //itmax = 3;
            use_xk = true;
            //verbose = true;
            //threshold = 1e-8;
            bcg(b);
            use_xk = false;

            for(std::map<int,int>::iterator it = glob_to_local.begin(); it != glob_to_local.end(); it++){

                if(it->first >= n_o && comm_map[it->second] == 1){
                    //if(inter_comm.rank() > my_bro) continue;

                    //cout << "Comm map : " << it->first << " " << comm_map[it->first] << endl;
                    vc.push_back(i);
                    vr.push_back(it->first - n_o);
                    vv.push_back(Xk(it->second,0));
                }
            }
        }
    }
    if(vv.size() == 0) {
        vc.push_back(0);
        vr.push_back(0);
        vv.push_back(0.0);
    }

    int *ii = new int[vv.size()];
    int *jj = new int[vv.size()];
    double *val = new double[vv.size()];
    for(int i = 0; i < vv.size(); i++){
        ii[i] = vr[i];
        jj[i] = vc[i];
        val[i] = vv[i];
    }
    shur = Coord_Mat_double(size_c, size_c, vv.size(), val, ii, jj);

    use_xk = false;

    return shur;
}		/* -----  end of function abcd::buildS  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::prodSv
 *  Description:  
 * =====================================================================================
 */
    MV_ColMat_double
abcd::prodSv ( MV_ColMat_double &V )
{
    MV_ColMat_double W(n, V.dim(1), 0);
    MV_ColMat_double b; // just for decoration!
    MV_ColMat_double R(V.dim(0), V.dim(1), 0);

    std::map<int,int>::iterator iti;
    for(int i = 0; i < V.dim(0); i++){
        iti = glob_to_local.find(n_o + i);

        if(iti != glob_to_local.end()){

            for(int j = 0; j < V.dim(1); j++){
                W(glob_to_local[n_o + i] , j) = V(i, j);
            }

        }
    }

    W = W - sumProject(0e0, b, 1e0, W);

    for(int i = 0; i < V.dim(0); i++){
        iti = glob_to_local.find(n_o + i);

        if(iti != glob_to_local.end()){

            for(int j = 0; j < V.dim(1); j++){
               R(i, j) = W(glob_to_local[n_o + i] , j);
            }

        }
    }
    return R;
}		/* -----  end of function abcd::prodSv  ----- */
