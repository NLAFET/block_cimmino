#include <abcd.h>
#include <iostream>
#include <fstream>


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::buildS
 *  Description:  Builds a sparse S
 * =====================================================================================
 */
    Coord_Mat_double
abcd::buildS (  )
{
    std::vector<int> cols;
    return buildS( cols );
}
    Coord_Mat_double
abcd::buildS ( std::vector<int> cols )
{
    std::vector<int> vc, vr;
    std::vector<double> vv;

    buildS(vr, vc, vv, cols);
    return Coord_Mat_double(size_c, size_c, vv.size(), &vv[0], &vr[0], &vc[0]);
}

void abcd::buildS(std::vector<int> &vr, std::vector<int> &vc, std::vector<double> &vv)
{
    std::vector<int> cols;
    buildS(vr, vc, vv, cols);
}

void abcd::buildS(std::vector<int> &vr, std::vector<int> &vc, std::vector<double> &vv, std::vector<int> &cols)
{
    Coord_Mat_double shur;
    //MV_ColMat_double S(size_c, size_c, 0);
    std::vector<int> sr, sc;
    std::vector<double> sv;

    std::vector<int> my_cols;

    if(cols.size() == 0) {
        for( int i = 0; i < size_c; i++){
            std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);
            if(iti!=glob_to_local.end()) my_cols.push_back(i);
        }
    } else {
        for( size_t i = 0; i < cols.size(); i++){
            std::map<int,int>::iterator iti = glob_to_local.find(n_o + cols[i]);
            if(iti!=glob_to_local.end()) my_cols.push_back(cols[i]);
        }
    }

    int maxcols = mpi::all_reduce(inter_comm, (int) my_cols.size(), mpi::maximum<int>());
    int mincols = mpi::all_reduce(inter_comm, (int) my_cols.size(), mpi::minimum<int>());
    int total = mpi::all_reduce(inter_comm, (int) my_cols.size(), std::plus<int>());
    IFMASTER clog << "Max number of cols is " << maxcols << " ,  min is " << mincols <<
       "and average is " << total/parallel_cg <<    endl;

#ifndef NO_MUMPS_ES
    mumps.keep[235 - 1] = icntl[Controls::exploit_sparcity];
    mumps.keep[261 - 1] = icntl[Controls::exploit_sparcity];
    mumps.keep[495 - 1] = icntl[Controls::exploit_sparcity];
    mumps.keep[497 - 1] = icntl[Controls::exploit_sparcity];
#endif
    if(dcntl[Controls::aug_type] == 0 || icntl[Controls::aug_iterative] == 2){

        std::vector<int>::iterator pos = my_cols.begin();
        std::vector<int>::iterator end_pos;

        vc.reserve(my_cols.size() * my_cols.size());
        vr.reserve(my_cols.size() * my_cols.size());
        vv.reserve(my_cols.size() * my_cols.size());

        int share = icntl[Controls::aug_blocking];

        while(pos != my_cols.end()){
            if(pos + share < my_cols.end()) end_pos = pos + share;
            else end_pos = my_cols.end();
            
            double perc = end_pos - my_cols.begin();
            perc /= my_cols.size();
            perc *= 100;


            std::vector<int> cur_cols;

            std::copy(pos, end_pos, std::back_inserter(cur_cols));


            //int mumps_share = share > 32 ? share : 16;
            int mumps_share = share;
            mumps.icntl[27 - 1] = mumps_share;

            // debug
            bool dense_build = false;
            if(dense_build){
                MV_ColMat_double sp = spSimpleProject(cur_cols);

                double *sptr = sp.ptr();
                int slda = sp.lda();
                int srows = sp.dim(0);

                for( size_t j = 0; j < cur_cols.size(); j++){
                    int c = cur_cols[j];
                    for( int i = c; i < srows; i++){
                        if(sptr[i + j * slda] != 0){
                            vr.push_back(i + 1);
                            vc.push_back(c + 1);
                            vv.push_back(sptr[i + j * slda]);
                        }
                    }
                }
            } else {
                spSimpleProject(cur_cols, vr, vc, vv);
            }
            pos = end_pos;
        }


    } else {

        for( int i = 0; i < size_c; i++){
            if(inter_comm.rank() == 0) clog << " Column " << i << " out of " << size_c << endl;

            icntl[Controls::block_size] = 1;
            Xk = MV_ColMat_double(n, 1, 0);
            MV_ColMat_double b(m, 1, 0); 

            std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);
            if(iti!=glob_to_local.end()){
                Xk(glob_to_local[n_o + i], 0) = 1;
            }

            use_xk = true;
            bcg(b);
            use_xk = false;

            for(std::map<int,int>::iterator it = glob_to_local.begin(); it != glob_to_local.end(); ++it){

                if(it->first >= n_o && comm_map[it->second] == 1){
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

    use_xk = false;

}       /* -----  end of function abcd::buildS  ----- */
