// Copyright Institut National Polytechnique de Toulouse (2014) 
// Contributor(s) :
// M. Zenadi <mzenadi@enseeiht.fr>
// D. Ruiz <ruiz@enseeiht.fr>
// R. Guivarch <guivarch@enseeiht.fr>

// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html"

// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 

// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 

// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.

#include <abcd.h>
#include <iostream>
#include <fstream>


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::buildS
 *  Description:  Builds a sparse S
 * =====================================================================================
 */
#ifdef WIP
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
#endif //WIP

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
    if(inter_comm.rank() == 0){
        LDEBUG << "Max number of cols is " << maxcols;
        LDEBUG << "Min number of cols is " << mincols;
        LDEBUG << "Avg number of cols is " << total/parallel_cg;
    }

#ifndef NO_MUMPS_ES
    mumps.keep[235 - 1] = icntl[Controls::exploit_sparcity];
    mumps.keep[261 - 1] = icntl[Controls::exploit_sparcity];
    mumps.keep[495 - 1] = icntl[Controls::exploit_sparcity];
    mumps.keep[497 - 1] = icntl[Controls::exploit_sparcity];
#endif

#ifdef WIP
    // If we need to fully augment the matrix, or at least to build part of it
    // this is needed only in ABCD direct and iterative
    if(dcntl[Controls::aug_type] == 0 || icntl[Controls::aug_iterative] == 2){
#endif

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

#ifdef WIP
            // debug
            bool dense_build = true;
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
#else
                spSimpleProject(cur_cols, vr, vc, vv);
#endif // WIP
            pos = end_pos;
        }


#ifdef WIP
    } else {
      // We build a smaller S from a filtered C, this requires an iterative
      // process. We repetedly compute m.n.s. using bcg() 

        for( int i = 0; i < size_c; i++){
            if(inter_comm.rank() == 0) LDEBUG << " Column " << i << " out of " << size_c;

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
#endif // WIP

    if(vv.size() == 0) {
        vc.push_back(0);
        vr.push_back(0);
        vv.push_back(0.0);
    }

    use_xk = false;

}       /* -----  end of function abcd::buildS  ----- */
