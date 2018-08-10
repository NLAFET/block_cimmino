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

/*!
 * \file sutils/buildS.cpp
 * \brief Implementation of the building of the sparse matrix S
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include <abcd.h>
#include <iostream>
#include <fstream>

#ifdef WIP
/*!
 *  \brief Build the matrix S directly or iteratively with no output
 *
 *  Build the matrix S directly or iteratively using sparse or dense RHS.
 *  The complete matrix will be computed but not output
 *
 *  \param vr: output array of rows of S
 *  \param vc: output array of columns of S
 *  \param vv: output array of values of S
 *
 */
    Coord_Mat_double
abcd::buildS (  )
{
    std::vector<int> cols;
    return buildS( cols );
}       /* -----  end of function abcd::buildS  ----- */

/*!
 *  \brief Build the matrix S directly or iteratively with no output
 *
 *  Build the matrix S directly or iteratively using sparse or dense RHS.
 *  A subpart can be selected but no output.
 *
 *  \param vr: output array of rows of S
 *  \param vc: output array of columns of S
 *  \param vv: output array of values of S
 *  \param cols: selected columns to compute only subpart of S
 *
 */
    Coord_Mat_double
abcd::buildS ( std::vector<int> cols )
{
    std::vector<int> vc, vr;
    std::vector<double> vv;

    buildS(vr, vc, vv, cols);
    return Coord_Mat_double(size_c, size_c, vv.size(), &vv[0], &vr[0], &vc[0]);
}       /* -----  end of function abcd::buildS  ----- */
#endif //WIP

/*!
 *  \brief Build the matrix S directly or iteratively
 *
 *  Build the matrix S directly or iteratively using sparse or dense RHS.
 *  The complete matrix will be computed.
 *
 *  \param vr: output array of rows of S
 *  \param vc: output array of columns of S
 *  \param vv: output array of values of S
 *
 */
void abcd::buildS(std::vector<int> &vr, std::vector<int> &vc, std::vector<double> &vv)
{
    std::vector<int> cols;
    buildS(vr, vc, vv, cols);
}       /* -----  end of function abcd::buildS  ----- */

/*!
 *  \brief Build the matrix S directly or iteratively
 *
 *  Build the matrix S directly or iteratively using sparse or dense RHS.
 *  A subpart can be selected.
 *
 *  \param vr: output array of rows of S
 *  \param vc: output array of columns of S
 *  \param vv: output array of values of S
 *  \param cols: selected columns to compute only subpart of S
 *
 */
void abcd::buildS(std::vector<int> &vr, std::vector<int> &vc, std::vector<double> &vv, std::vector<int> &cols)
{
    Coord_Mat_double shur;
    //MV_ColMat_double S(size_c, size_c, 0);

    // arrays of rows/columns/values of the matrix S
    std::vector<int> sr, sc;
    std::vector<double> sv;

    // find local columns
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

    // Display maximum/minimum/total number of columns
    int maxcols = mpi::all_reduce(inter_comm, (int) my_cols.size(), mpi::maximum<int>());
    int mincols = mpi::all_reduce(inter_comm, (int) my_cols.size(), mpi::minimum<int>());
    int total = mpi::all_reduce(inter_comm, (int) my_cols.size(), std::plus<int>());
    if(inter_comm.rank() == 0){
        LDEBUG << "Max number of cols is " << maxcols;
        LDEBUG << "Min number of cols is " << mincols;
        LDEBUG << "Avg number of cols is " << total/parallel_cg;
    }

#ifdef WIP

    // Mumps exploi sparsity ?????
#ifndef NO_MUMPS_ES
    mumps.keep[235 - 1] = icntl[Controls::exploit_sparcity];
    mumps.keep[261 - 1] = icntl[Controls::exploit_sparcity];
    mumps.keep[495 - 1] = icntl[Controls::exploit_sparcity];
    mumps.keep[497 - 1] = icntl[Controls::exploit_sparcity];
#endif //no mumps exploits parsity

    // If we need to fully augment the matrix, or at least to build part of it
    // this is needed only in ABCD direct and iterative
    if(icntl[Controls::aug_type] == 0 || icntl[Controls::aug_iterative] == 2){
#endif //wip
        std::vector<int>::iterator pos = my_cols.begin();
        std::vector<int>::iterator end_pos;

        // Allocate memory as if local part of S is dense
        vc.reserve(my_cols.size() * my_cols.size());
        vr.reserve(my_cols.size() * my_cols.size());
        vv.reserve(my_cols.size() * my_cols.size());

        // number of columns to compute at the same time (chunk)
        int share = icntl[Controls::aug_blocking];

        while(pos != my_cols.end()){
            // if the next chunk does not arrive at the end
            if(pos + share < my_cols.end()) end_pos = pos + share;
            else end_pos = my_cols.end();

            // percentage of columns computed
            double perc = end_pos - my_cols.begin();
            perc /= my_cols.size();
            perc *= 100;

            // current local columns
            std::vector<int> cur_cols;
            std::copy(pos, end_pos, std::back_inserter(cur_cols));

            mumps.icntl[27 - 1] = share; // blocking size for multiple RHS

#ifdef WIP
            // debug
            bool dense_build = false;
            // use simple_sumproject.cpp
            if(dense_build){
                MV_ColMat_double sp = spSimpleProject(cur_cols);

                double *sptr = sp.ptr(); // pointer to S part
                int slda = sp.lda(); // columns of S part
                int srows = sp.dim(0); // rows of S part

                // Save S in global arrays
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
            // use sparse_simple_sumproject.cpp
            } else {
#endif // WIP
                spSimpleProject(cur_cols, vr, vc, vv);
#ifdef WIP
            }
#endif // WIP
            pos = end_pos;
        }


#ifdef WIP
    } else {
      // We build a smaller S from a filtered C, this requires an iterative
      // process. We repetedly compute m.n.s. using bcg()
        for( int i = 0; i < size_c; i++){
            if(inter_comm.rank() == 0) LDEBUG << " Column " << i << " out of " << size_c;

            // Initialize Xk (0) and b (0) for BCG
            icntl[Controls::block_size] = 1;
            Xk = MV_ColMat_double(n, 1, 0);
            MV_ColMat_double b(m, 1, 0);

            // Set local columns of interest to 1
            std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);
            if(iti!=glob_to_local.end()){
                Xk(glob_to_local[n_o + i], 0) = 1;
            }

            // BCG with starting point Xk
            use_xk = true;
            bcg(b);
            use_xk = false;

            // Save computed S
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

    // If S has size 0, then it is 0
    if(vv.size() == 0) {
        vc.push_back(0);
        vr.push_back(0);
        vv.push_back(0.0);
    }

    use_xk = false;
}       /* -----  end of function abcd::buildS  ----- */
