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
 * \file comms.cpp
 * \brief Implementation of the creation of intercommunicators and distribution of RHS
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include <abcd.h>
#include <vect_utils.h>

/*!
 *  \brief Assigns each mpi-process to its category : CG-master or MUMPS-Slave
 *
 *  Decides which processes are Masters and create communicators between masters
 *  and between slaves. The masters are decided:
 *    - choice 1: zig-zag distribution, 1master-1node
 *    - choice -1: first processes are masters (Momo's implementation)
 *
 */
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

    // Broadcast choice for Master/Slave distribution
    mpi::broadcast(comm, icntl[Controls::master_def], 0);
    mpi::broadcast(comm, icntl[Controls::slave_def], 0);

    // Get communication map information on root MPI
    mpi::environment env;
    std::string node = env.processor_name();
    int cpu = sched_getcpu();
    int root_node;
    std::vector<int> node_count; // number of MPI per node
    std::vector<std::pair<int, std::string>> pair_vect; // vector of pairs
    std::pair<int, std::string> p (comm.rank(), node); // pair mpi-node
    mpi::gather(comm, p, pair_vect, 0);
    if (comm.rank() == 0) {
        int count=0; // number of nodes
        std::map<std::string, int> nodes; // index of nodes
        for(int iii=0; iii<pair_vect.size(); ++iii) {
            // if we meet the node for the first time
            if (nodes.find(pair_vect[iii].second) == nodes.end()) {
                nodes[pair_vect[iii].second] = count;
                node_count.push_back(0);
                std::vector<int> v;
                node_map_slaves.push_back(v);
                ++count;
            }
            // increase number of MPI in this node
            ++node_count[nodes[pair_vect[iii].second]];
            if (pair_vect[iii].first)
                node_map_slaves[nodes[pair_vect[iii].second]].push_back(iii);
            // save the node index for this MPI
            mpi_map.push_back(nodes[pair_vect[iii].second]);
        }
        root_node=nodes[node];
        // now the vector nodes is useless
    }

    int choice=icntl[Controls::master_def];
    // choice 1: zig-zag 1master/node assignation
    if (choice == 1) {
        // This choice will place 1 master on each node in descending order or node size
        // then if masters remain will continue adding them in zig-zag to the nodes
        if (comm.rank() == 0) {
	    // Sort nodes per count via tmp_ord to ord_nodes
            // The node of the root MPI must always be first
            std::vector<std::pair<int, int>> tmp_ord; // tmp vector to sort nodes by #MPI
            for (int iii=0; iii<node_count.size(); ++iii) {
                if (iii != root_node) {
                    std::pair<int, int> p (node_count[iii], iii);
                    tmp_ord.push_back(p);
                }
            }
            std::sort(tmp_ord.begin(), tmp_ord.end(),
                pair_comparison<int, int,
                    position_in_pair::first,
                    comparison_direction::descending>);
            std::vector<int> ord_nodes; // node ordered index vs. original index
            ord_nodes.push_back(root_node);
            for (int iii=0; iii<tmp_ord.size(); ++iii) {
                ord_nodes.push_back(tmp_ord[iii].second);
            }

            // Assign masters from biggest node to smallest then from smallest to biggest etc. (zig-zag)
            instance_type_vect.assign(comm.size(), 1);
            if (node_map_slaves.size() > 1) {
                int current_idx=1;
                int direction=1;
                masters_node.push_back(0);
                instance_type_vect[0] = 0;
                for (int iii=1; iii<parallel_cg; ++iii) {
                    int current_node=ord_nodes[current_idx];
                    // check if remaining MPI in node and assign as master
		    if (node_map_slaves[current_node].size() > 0) {
                        masters_node.push_back(current_node);
                        instance_type_vect[node_map_slaves[current_node].back()] = 0;
                        node_map_slaves[current_node].pop_back();
                    } else --iii; //do not assign, only update the zig-zag direction
                    // update direction of node exploration
                    if (!current_idx) {
                        ++direction;
                    } else if (!((current_idx+1)%node_map_slaves.size())) {
                        --direction;
                    }
                    current_idx+=direction;
                }
            // If only 1 node, all masters in first node
            } else {
                masters_node.push_back(0);
                instance_type_vect[0] = 0;
                for (int iii=1; iii<parallel_cg; ++iii) {
                    masters_node.push_back(0);
                    instance_type_vect[node_map_slaves[0].back()] = 0;
                    node_map_slaves[0].pop_back();
                }
            }
        }

        mpi::broadcast(comm, instance_type_vect, 0);
        instance_type = instance_type_vect[comm.rank()];
    // choice -1: previous implementation
    } else {
        instance_type = comm.rank() < parallel_cg ? 0 : 1;
        if (comm.rank() == 0) {
            for(int i=0; i<parallel_cg; ++i) {
                masters_node.push_back(mpi_map[i]);
            }
        }
    }

    // Define communicators between masters and between slaves
    inter_comm = comm.split(instance_type == 0);

    if (comm.rank() == 0)
        LINFO << "There are " << inter_comm.size()
              << " master processes and "
              << comm.size() - inter_comm.size()
              << " workers";
}               /* -----  end of function abcd::createInterCommunicators  ----- */

/*!
 *  \brief Distribute the local RHS to masters
 *
 *  Distribute local parts of the RHS corresponding to each partition for each
 *  of the block_size vectors.
 *
 */
void abcd::distributeRhs()
{
    // Broadcast if artificial RHS and Block Size (=max(block Size, mrhs))
    mpi::broadcast(inter_comm, use_xf, 0);
    if(icntl[Controls::block_size] < nrhs) icntl[Controls::block_size] = nrhs;
    mpi::broadcast(inter_comm, icntl[Controls::block_size], 0);

    if(comm.rank() == 0) {
        /* Build the RHS */
        // size of RHS = accumulated size of partitions
        int r = 0;
        for(int i = 0; i < nb_local_parts; i++){
            r += partitions[i].dim(0);
        }

        // if artificial RHS used, computed its Inf-norm
	if(Xf.dim(0) != 0) {
	    nrmXf = 0;
	    for(int i = 0; i < A.dim(1); i++) {
  	        if(abs(Xf(i,0)) > nrmXf) nrmXf = abs(Xf(i,0));
	    }
	}
        B = MV_ColMat_double(m_l, icntl[Controls::block_size], 0);
        for(int j = 0; j < nrhs; j++){
            for(int i = 0; i < m_l; i++) {
                B(i, j) = rhs[i + j*m_l] * drow_[i];
            }
        }

        // Check RHS is not zero and broadcast
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

        // If the Block Size is bigger than the nrhs, generate additional random vectors with values in [1,2]
        if(icntl[Controls::block_size] > nrhs) {
            // Create additional RHS block BR
            double *rdata = new double[n_l * (icntl[Controls::block_size] - nrhs)];
            srand(n_l);
            for(int i=0; i< n_l*(icntl[Controls::block_size]-nrhs); i++){
                rdata[i] = (double)((rand())%100+1)/99.9 + 1;
            }
            MV_ColMat_double BR(rdata, n_l, icntl[Controls::block_size] - nrhs, MV_Matrix_::ref);

            // Add this part multiplied by A to B (B=[B RR]=[B A*BR])
            MV_ColMat_double RR = smv(A, BR);
            B(MV_VecIndex(0,B.dim(0)-1),MV_VecIndex(nrhs,icntl[Controls::block_size]-1)) =
                RR(MV_VecIndex(0,B.dim(0)-1), MV_VecIndex(0, icntl[Controls::block_size]-nrhs - 1));
            delete[] rdata;
        }

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

        // For each master, send the corresponding part of the RHS (bi of partitions in each vector)
        for(int k = 1; k < parallel_cg; k++) {
            for(int j = 0; j < icntl[Controls::block_size]; j++) {
                for(size_t i = 0; i < partitionsSets[k].size(); i++){
                    int p = partitionsSets[k][i];

		    // Not tested for blocksize >1
		    double *tmp = (double*) malloc(sizeof(double)*row_indices[p].size()*icntl[Controls::block_size]);
		    for(int zz=0;  zz< row_indices[p].size(); zz++ ){
			tmp[zz] = b_ptr[row_indices[p][zz] + j * m_l ];
		    }
                    inter_comm.send(k, 18, tmp, row_indices[p].size());
                }
            }
        }

        // If no artificial RHS, restrict the local RHS to only the local partitions
        if(!use_xf){
            MV_ColMat_double BB(B);

            B = MV_ColMat_double(m, icntl[Controls::block_size], 0);
            int pos = 0;
            for(size_t i = 0; i < partitionsSets[0].size(); i++){
                int p = partitionsSets[0][i];
                for(int j = 0; j < row_indices[p].size() ; j++){
                    for(int k = 0; k < icntl[Controls::block_size]; k++){
                        B(pos, k) = BB( row_indices[p][j], k);
		    }
		    pos++;
		}
            }
        }
    } else {
        // Check if good RHS from root
        int good_rhs;
        mpi::broadcast(inter_comm, good_rhs, 0);
        if (good_rhs != 0) {
            info[Controls::status] = good_rhs;
            stringstream err_msg;
            err_msg << "On process [" << comm.rank() << "], leaving due to an error on the master";
            throw std::runtime_error(err_msg.str());
        }

        // Receive temp solution
        Xf = MV_ColMat_double(n_o, 1, 0);
        double *xf_ptr = Xf.ptr();
        inter_comm.recv(0, 171, xf_ptr, n_o);

        /* Receive RHS */
        //Initialize RHS with zeros
        inter_comm.recv(0, 17, nrhs);
        int size_rhs = m*icntl[Controls::block_size];
        rhs = new double[size_rhs];
        for(int i = 0; i < m * icntl[Controls::block_size]; i++) rhs[i] = 0;

        // Receive and save local parts of the RHS (for each partition)
        B = MV_ColMat_double(m, icntl[Controls::block_size], 0);
        for(int j = 0; j < icntl[Controls::block_size]; j++) {
            int p = 0;
            // receive current vector bi for each partition
            for(int k = 0; k < nb_local_parts; k++){
                inter_comm.recv(0, 18, rhs + p + j * m, partitions[k].dim(0));
                p+= partitions[k].dim(0);
            }

            // save current RHS vector
            VECTOR_double t(rhs+j*m, m);
            B.setCol(t, j);
        }

        delete[] rhs;
    }
    // Broadcast max iterations/threshold for BCG/threshold for filtering augmentation
    mpi::broadcast(inter_comm, icntl[Controls::itmax], 0);
    mpi::broadcast(inter_comm, dcntl[Controls::threshold], 0);
#ifdef WIP
    mpi::broadcast(inter_comm, dcntl[Controls::aug_filter], 0);
#endif //WIP
}               /* -----  end of function abcd::distributeRhs  ----- */
