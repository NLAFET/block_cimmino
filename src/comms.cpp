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
#include <vect_utils.h>

/// Assignes each mpi-process to its category : CG-master or MUMPS-Slave
void abcd::createInterCommunicators()
{
    mpi::group grp = comm.group();

    mpi::communicator cm;
    mpi::group gp;

    mpi::broadcast(comm, parallel_cg, 0);

    if(parallel_cg > comm.size()){
        info[Controls::status] = -8;
        throw std::runtime_error("The number of masters is larger than the number of MPI-processes");
    }

    instance_type = comm.rank() < parallel_cg ? 0 : 1;

    inter_comm = comm.split(instance_type == 0);

    if (comm.rank() == 0)
        LINFO << "There are " << inter_comm.size()
              << " master processes and "
              << comm.size() - inter_comm.size()
              << " workers";
}

void abcd::distributeRhs()
{
    mpi::broadcast(inter_comm, use_xf, 0);

    if(icntl[Controls::block_size] < nrhs) icntl[Controls::block_size] = nrhs;
    mpi::broadcast(inter_comm, icntl[Controls::block_size], 0);

    if(comm.rank() == 0) {

        int r_pos = 0;
        // Build my part of the RHS
        //int r = std::accumulate(partitions.begin(), partitions.end(), 0, sum_rows);
        int r = 0;
        for(int i = 0; i < nb_local_parts; i++){
            r += partitions[i].dim(0);
        }
        
        if(rhs == nullptr){
            rhs = new double[n_l * nrhs];

            srand(10); 
            B = MV_ColMat_double(m_l, icntl[Controls::block_size]);

            nrmXf = 0;
            Xf = MV_ColMat_double(A.dim(1), nrhs);
            for(int j = 0; j < nrhs; j++){
                for(int i = 0; i < A.dim(1); i++){
                    rhs[i + j * n_l] = (double)((rand())%100+1)/99.0;
                }
            }

            for(int j = 0; j < nrhs; j++){
                VECTOR_double xf_col(A.dim(1));
                for(int i = 0; i < A.dim(1); i++) {
                    xf_col[i] = rhs[i + j * A.dim(1)];
                }
                Xf.setCol(xf_col, j);
            }

            MV_ColMat_double BB = smv(A, Xf);

            for(int j = 0; j < nrhs; j++){
                double unscaled; 
                for(int i = 0; i < A.dim(1); i++) {
                    unscaled = rhs[i + j * A.dim(1)] * dcol_[i];
                    if(abs(unscaled) > nrmXf) nrmXf = abs(unscaled);
                    Xf(i, j) = unscaled;
                }
            }

            for(int j = 0; j < nrhs; j++){
                //VECTOR_double t(rhs+j*m_l, m_l);
                //B.setCol(t, j);
                //B.push_back(t);
                B(MV_VecIndex(0, m_l-1), MV_VecIndex(0,nrhs-1)) = BB;
            }
        } else {
            B = MV_ColMat_double(m_l, icntl[Controls::block_size], 0);
            if(row_perm.size() != 0){
                for(int j = 0; j < nrhs; j++){
                    for(int i = 0; i < m_l; i++) {
                        B(i, j) = rhs[row_perm[i] + j*m_l];
                    }
                }

            } else {
                for(int j = 0; j < nrhs; j++){
                    for(int i = 0; i < m_l; i++) {
                        B(i, j) = rhs[i + j*m_l];
                    }
                }
            }

            diagScaleRhs(B);

        }

        int good_rhs = 0;
        if (infNorm(B) == 0) {
            good_rhs = -9;
            mpi::broadcast(inter_comm, good_rhs, 0);
            stringstream err_msg;
            err_msg << "On process [" << comm.rank() << "], the given right-hand side is zero";
            info[Controls::status] = good_rhs;
            throw std::runtime_error(err_msg.str());
        }
        
        mpi::broadcast(inter_comm, good_rhs, 0);

        if(icntl[Controls::block_size] > nrhs) {
            double *rdata = new double[n_l * (icntl[Controls::block_size] - nrhs)];

            srand(n_l); 
            for(int i=0; i< n_l*(icntl[Controls::block_size]-nrhs); i++){ 
                rdata[i] = (double)((rand())%10)/99.9 + 1;
                //rdata[i] = i+1;
            }
            MV_ColMat_double BR(rdata, n_l, icntl[Controls::block_size] - nrhs, MV_Matrix_::ref);
            MV_ColMat_double RR = smv(A, BR);

            B(MV_VecIndex(0,B.dim(0)-1),MV_VecIndex(nrhs,icntl[Controls::block_size]-1)) = 
                RR(MV_VecIndex(0,B.dim(0)-1), MV_VecIndex(0, icntl[Controls::block_size]-nrhs - 1));
            delete[] rdata;
        }

        r_pos += r;

        double *b_ptr = B.ptr();

        // temp solution :
        // send Xf to everybody
        double *xf_ptr = Xf.ptr();
        // for other masters except me!
        for(int k = 1; k < parallel_cg; k++) {
            inter_comm.send(k, 171, xf_ptr, Xf.dim(0));
        }

        // for other masters except me!
        for(int k = 1; k < parallel_cg; k++) {
            // get the partitions that will be sent to the master
            inter_comm.send(k, 17, nrhs);
        }

        for(int k = 1; k < parallel_cg; k++) {
            for(int j = 0; j < icntl[Controls::block_size]; j++) {
                for(size_t i = 0; i < partitionsSets[k].size(); i++){
                    int p = partitionsSets[k][i];
                    inter_comm.send(k, 18, b_ptr + strow[p] + j * m_l, nbrows[p]);
                }
            }
        }

        if(!use_xf){
            MV_ColMat_double BB(B);

            B = MV_ColMat_double(m, icntl[Controls::block_size], 0);
            int pos = 0;
            for(size_t i = 0; i < partitionsSets[0].size(); i++){
                int p = partitionsSets[0][i];
                for(int j = 0; j < nbrows[p]; j++){
                    for(int k = 0; k < icntl[Controls::block_size]; k++){
                        B(pos, k) = BB(strow[p] + j, k);
                    }
                    pos++;
                }
            }
        }
    } else {
        int good_rhs;
        mpi::broadcast(inter_comm, good_rhs, 0);
        if (good_rhs != 0) {
            info[Controls::status] = good_rhs;
            stringstream err_msg;
            err_msg << "On process [" << comm.rank() << "], leaving due to an error on the master";
            throw std::runtime_error(err_msg.str());
        }

        Xf = MV_ColMat_double(n_o, 1, 0);
        double *xf_ptr = Xf.ptr();
        inter_comm.recv(0, 171, xf_ptr, n_o);

        inter_comm.recv(0, 17, nrhs);

        int size_rhs = m*icntl[Controls::block_size];
        rhs = new double[size_rhs];
        for(int i = 0; i < m * icntl[Controls::block_size]; i++) rhs[i] = 0;

        B = MV_ColMat_double(m, icntl[Controls::block_size], 0);
        for(int j = 0; j < icntl[Controls::block_size]; j++) {
            int p = 0;
            for(int k = 0; k < nb_local_parts; k++){
                inter_comm.recv(0, 18, rhs + p + j * m, partitions[k].dim(0));
                p+= partitions[k].dim(0);
            }

            VECTOR_double t(rhs+j*m, m);
            B.setCol(t, j);
        }

        delete[] rhs;
    }
    // and distribute max iterations
    mpi::broadcast(inter_comm, icntl[Controls::itmax], 0);
    mpi::broadcast(inter_comm, dcntl[Controls::threshold], 0);
#ifdef WIP
    mpi::broadcast(inter_comm, dcntl[Controls::aug_filter], 0);
#endif //WIP

    // A = CompRow_Mat_double();
}
