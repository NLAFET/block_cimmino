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

void abcd::distributeData()
{
    if(comm.rank() == 0) {
        std::vector<int> nnz_parts;
        std::vector<int> m_parts;
        std::vector<int> groups;

        for(int k = 0; k < icntl[Controls::nbparts]; k++) {
            nnz_parts.push_back(parts[k].NumNonzeros());
            // Criteria for partition distribution
            m_parts.push_back(parts[k].dim(0));
        }

        abcd::partitionWeights(partitionsSets, m_parts, parallel_cg);

        for(int i = 1; i < parallel_cg ; i++) {
            inter_comm.send(i, 0, partitionsSets[i]);

            for(unsigned int k = 0; k < partitionsSets[i].size(); k++){
                int j = partitionsSets[i][k];

                std::vector<int> dims;
                dims.push_back(parts[j].dim(0));
                dims.push_back(parts[j].dim(1));

		// send the nnz
                inter_comm.send(i, 1, parts[j].NumNonzeros());
		// send the dimensions
                inter_comm.send(i, 2, dims);
		// send the origin number of columns of the matrix
                inter_comm.send(i, 21, n);

		// send the partitions data
                inter_comm.send(i, 3, parts[j].colind_ptr(), parts[j].NumNonzeros());
                inter_comm.send(i, 4, parts[j].rowptr_ptr(), dims[0] + 1);
                inter_comm.send(i, 5, parts[j].val_ptr(), parts[j].NumNonzeros());

		// send the column index data
                inter_comm.send(i, 6, column_index[j]);

		// if we augment the matrix, send the indices where C starts
                if(icntl[Controls::aug_type] > 0) inter_comm.send(i, 61, stC[j]);
            }

        }
        LINFO << "Sent partitions";

        m_l = m;
        n_l = n;

        if(parallel_cg != 1){
            std::vector<std::vector<int> > columnIndices;
            std::vector<int> stcs;

	    // duplicate my partitions data so that we clear other 
	    // processes data from the master's memory
            for(unsigned int i = 0; i < partitionsSets[0].size(); i++){
                int j = partitionsSets[0][i];
                partitions.push_back(parts[j]);
                m += parts[j].dim(0);
                nz += parts[j].NumNonzeros();

                columnIndices.push_back(column_index[j]);
                if(icntl[Controls::aug_type] != 0) stcs.push_back(stC[j]);
            }
            parts.clear();
            column_index.clear();
            nb_local_parts = partitions.size();
	    
            if(icntl[Controls::aug_type] != 0) stC.clear();
            for(int i = 0; i < nb_local_parts; i++){
                //partitions[i] = tp[i];
                column_index.push_back(columnIndices[i]);
                if(icntl[Controls::aug_type] != 0) stC.push_back(stcs[i]);
            }
        } else {
            for(unsigned int i = 0; i < parts.size(); i++){
                partitions.push_back(parts[i]);
//	        partitionsSets[0][i]=i;
            }
            parts.clear();
            nb_local_parts = partitions.size();
        }
        m = 0;
        nz = 0;

        for(int i = 0; i < nb_local_parts; i++){
            m += partitions[i].dim(0);
            nz += partitions[i].NumNonzeros();
        }


    } else {
        std::vector<int> se;
        inter_comm.recv(0, 0, se);
        int sm = 0, snz = 0;
        for(unsigned int i = 0; i < se.size(); i++) {


            // THe partition data
            std::vector<int> dims;
            int l_nz, l_m, l_n;
            inter_comm.recv(0, 1, l_nz);
            inter_comm.recv(0, 2, dims);
            inter_comm.recv(0, 21, n);
            l_m = dims[0];
            l_n = dims[1];

            int *l_jcn = new int[l_nz];
            int *l_irst = new int[l_m + 1];
            double *l_v = new double[l_nz];

            inter_comm.recv(0, 3, l_jcn, l_nz);
            inter_comm.recv(0, 4, l_irst, l_m + 1);
            inter_comm.recv(0, 5, l_v, l_nz);

            // Continue the communications and let the initialization of the matrix at the end
            std::vector<int> ci;
            inter_comm.recv(0, 6, ci);
            column_index.push_back(ci);

            int stc;
            if(icntl[Controls::aug_type] > 0){
                inter_comm.recv(0, 61, stc);
                stC.push_back(stc);
            }


            // Create the matrix and push it in!
            CompRow_Mat_double mat(l_m, l_n, l_nz, l_v, l_irst, l_jcn);
            partitions.push_back(mat);

            sm += l_m;
            snz += l_nz;
            delete[] l_jcn;
            delete[] l_irst; 
            delete[] l_v;
        }
        nb_local_parts = partitions.size();
        LDEBUG3 << "Process " << inter_comm.rank() << " received " << nb_local_parts << " partitions";

        // Set the number of rows and nnz handled by this CG Instance
        m = sm;
        nz = snz;
    }
    mpi::broadcast(inter_comm, m_l, 0);
    mpi::broadcast(inter_comm, n_l, 0);

    // compute the matrix inf-norm in a distributed maner
    nrmMtx = 0;
    double nrmP = 0;
    for(int i = 0; i < partitions.size(); ++i) {
      int *rp = partitions[i].rowptr_ptr();
      int *cp = partitions[i].colind_ptr();
      double *vp = partitions[i].val_ptr();
      
      for(int r = 0; r < partitions[i].dim(0); r++) {
          double rsum = 0;
          for (int c = rp[r]; c < rp[r+1]; ++c){
              rsum += abs(vp[c]);
          }
          if(nrmP < rsum) nrmP = rsum;
      }
    }
    mpi::all_reduce(inter_comm, &nrmP, 1, &nrmMtx, mpi::maximum<double>());
}

void abcd::createInterconnections()
{
    if(comm.rank() == 0) LINFO << "Creating interconnections between processes";
    
    // Link between the current partition and the global array
    // is used only in ABCD
    if (icntl[Controls::aug_type] != 0) {
        for(int k = 0; k < nb_local_parts; k++) {
            std::map<int, int> gt;
            std::map<int, int> ptg;
            for(unsigned int j = 0; j < column_index[k].size(); j++) {
                gt[column_index[k][j]] = j;
                ptg[j] = column_index[k][j];
            }
            part_to_glob.push_back(ptg);
            glob_to_part.push_back(gt);
        }
    }

    // we need the merge of column indices 
    std::vector<int> merge_index = mergeSortedVectors(column_index);

    // for ABCD, we need a global to local indices so that we can
    // identify which column in C is linked to 
    for(unsigned int j = 0; j < merge_index.size(); j++) {
        glob_to_local[merge_index[j]] = j;
        glob_to_local_ind.push_back(merge_index[j]);
    }
    if (icntl[Controls::aug_type] != 0) {
        // defines the starting point of C in the local columns
        st_c_part_it = glob_to_local_ind.end();

        while(*(st_c_part_it - 1) >= n_o) st_c_part_it--;

        st_c_part = st_c_part_it - glob_to_local_ind.begin();
    }

    // for each partition find a local column index for the previous merge
    local_column_index = std::vector<std::vector<int> >(nb_local_parts);

    for(int p = 0; p < nb_local_parts; p++) {
        unsigned int i = 0, j = 0;
        while (i < column_index[p].size() && j < merge_index.size()) {

            if(column_index[p][i] == merge_index[j]) {
                local_column_index[p].push_back(j);
                i++;
            }

            j++;
        }
    }


    // All-to-all exchange of column indices
    std::map<int, std::vector<int> > their_cols;
    std::vector<mpi::request> reqs_c;
    for (int i = 0; i < parallel_cg; i++) {
        if (i == inter_comm.rank())
            continue;
        reqs_c.push_back(inter_comm.irecv(i, 41, their_cols[i]));
        reqs_c.push_back(inter_comm.isend(i, 41, merge_index));
    }

    mpi::wait_all(reqs_c.begin(), reqs_c.end());

    for(std::map<int, std::vector<int> >::iterator it = their_cols.begin();
            it != their_cols.end(); ++it) {

        col_interconnections[it->first] =
            getIntersectionIndices(
                merge_index, their_cols[it->first]
            ).first;

    }

    // Create a Communication map for ddot
    comm_map.assign(n, 1);
    for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
            it != col_interconnections.end(); ++it) {
        // if I share data with it->first and I'm after him, let him compute!
        if(inter_comm.rank() > it->first) {
            for(std::vector<int>::iterator i = it->second.begin(); i != it->second.end(); ++i) {
                if(comm_map[*i] == 1) comm_map[*i] = -1;
            }
        }
    }

    if (inter_comm.rank() == 0) 
        LINFO << "Interconnections created";
    
}
