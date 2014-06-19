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

#include<abcd.h>
#include<mumps.h>
using namespace  boost::lambda;

#ifdef WIP
MV_ColMat_double abcd::spSimpleProject(std::vector<int> mycols)
{
    bool dense_rhs = (icntl[Controls::aug_dense] == 1);
    //dense_rhs = true;

    int s = mycols.size();
    std::vector<int> rr;
    std::vector<double> rv;

    CompRow_Mat_double *r = new CompRow_Mat_double[nb_local_parts];
    std::vector<std::map<int,int> > loc_cols(nb_local_parts);

    int nzr_estim = 0;

    for(int k = 0; k < nb_local_parts; k++) {

        CompRow_Mat_double Y;

        VECTOR_int yr(mycols.size(), 0);
        VECTOR_int yc(mycols.size(), 0);
        VECTOR_double yv(mycols.size(), 0);
        int c;

        int ct = 0;
        for(size_t i = 0; i < mycols.size(); i++){
            c = mycols[i];
            if(glob_to_part[k].find(n_o + c) != glob_to_part[k].end()){
                yr[ct] = glob_to_part[k][n_o + c];
                yc[ct] = i;
                yv[ct] = 1;

                ct++;
                loc_cols[k][i] = 1;
            }
        }

        Coord_Mat_double Yt(partitions[k].dim(1), s, ct, yv.ptr(), yr.ptr(), yc.ptr());

        Y = CompRow_Mat_double(Yt);

        r[k] = spmm(partitions[k], Y) ;
        nzr_estim += r[k].NumNonzeros();
    }


    if(dense_rhs){
        // Build the mumps rhs
        mumps.rhs = new double[mumps.n * s];
        for(int i = 0; i < mumps.n * s; i++) mumps.rhs[i] = 0;
        MV_ColMat_double mumps_rhs(mumps.rhs, mumps.n, s, MV_Matrix_::ref);

        // dense mumps rhs
        mumps.icntl[20 - 1] = 0;

        int pos = 0;
        for(int k = 0; k < nb_local_parts; k++) {
            for(int r_p = 0; r_p < s; r_p++) {
                for(int i = pos; i < pos + partitions[k].dim(1); i++) {
                    mumps_rhs(i, r_p) = 0;
                }
                int j = 0;
                for(int i = pos + partitions[k].dim(1); i < pos + partitions[k].dim(1) + partitions[k].dim(0); i++) {
                    mumps_rhs(i, r_p) = r[k](j, r_p);
                    j++;
                }
            }
            pos += partitions[k].dim(1) + partitions[k].dim(0);
        }
    } else {
        // sparse mumps rhs
        mumps.setIcntl(20, 1);

        // Build the mumps rhs
        mumps.rhs = new double[mumps.n * s];
        for(int i = 0; i < mumps.n * s; i++) mumps.rhs[i] = 0;
        MV_ColMat_double mumps_rhs(mumps.rhs, mumps.n, s, MV_Matrix_::ref);

        {
            mumps.irhs_ptr      = new int[s + 1];
            rr.reserve(nzr_estim);
            rv.reserve(nzr_estim);

            int cnz = 1;


            for(int r_p = 0; r_p < s; r_p++) {
                mumps.irhs_ptr[r_p] = cnz;

                int pos = 0;
                for(int k = 0; k < nb_local_parts; k++) {
                    CompRow_Mat_double rtt = r[k];
                    int _dim1 = partitions[k].dim(1);
                    int _dim0 = partitions[k].dim(0);

                    int j = 0;
                    for(int i = pos + _dim1; i < pos + _dim1 + _dim0; i++) {
                        if(rtt(j, r_p) != 0){
                            rr.push_back(i+1);
                            rv.push_back(rtt(j, r_p));
                            cnz++;
                        }
                        j++;
                    }
                    pos += _dim1 + _dim0;
                }
            }
            mumps.irhs_ptr[s] = cnz;

            mumps.nz_rhs        = cnz - 1;

            mumps.rhs_sparse    = &rv[0];
            mumps.irhs_sparse   = &rr[0];
        }

    }

    int job = 1;
    mpi::broadcast(intra_comm, job, 0);

    mumps.nrhs          = s;
    mumps.lrhs          = mumps.n;
    mumps.job           = 3;

    dmumps_c(&mumps);

    MV_ColMat_double Delta;

    Delta = MV_ColMat_double(size_c, s, 0);

    int dlda = Delta.lda();
    double *dpt = Delta.ptr();

    int x_pos = 0;

    for(int k = 0; k < nb_local_parts; k++) {
        int start_c;

        // find the begining of the C part, if there is no C, set it to the end of the current part
        if (stC[k] != -1) {
            start_c = glob_to_part[k][stC[k]];
        } else {
            start_c = column_index[k].size();
        }

        // move the pointer to the end of the current part (excluding C part)
        x_pos += start_c;

        // the actual sumproject, only on the C part, as it's the needed sum
        for(size_t i = start_c; i < column_index[k].size(); i++){
            int ci = column_index[k][i] - n_o;

            for(int j = 0; j < s; j++)
                dpt[ci + j * dlda] -= mumps.rhs[x_pos + j * mumps.n];

            x_pos++;
        }

        for(int j = 0; j < s; j++) {
            if(loc_cols[k][j]){
                int c = mycols[j];

                // as we have two partitions, add a half on each to obtain identity
                dpt[c + j * dlda] += 0.5;
            }
        }

        // move to the next partition
        x_pos += partitions[k].dim(0);
    }

    // disable sparse mumps rhs
    mumps.icntl[20 - 1] = 0;

    delete[] mumps.rhs;
    delete[] r;

    return Delta;
}
#endif // WIP
