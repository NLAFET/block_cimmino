#include <abcd.h>

/// Assignes each mpi-process to its category : CG-master or MUMPS-Slave
void abcd::createInterComm()
{
    mpi::communicator world;
    mpi::group grp = world.group();

    mpi::communicator cm;
    mpi::group gp;

    mpi::broadcast(world, parallel_cg, 0);

    if(parallel_cg > world.size()) throw - 14;

    inter_comm = world.split(world.rank() < parallel_cg);

    if(world.rank() < parallel_cg)
        instance_type = 0;
    else
        instance_type = 1;

}

/// Distribute the partitions over CG processes
void abcd::distributePartitions()
{
    mpi::communicator world;

    if(world.rank() == 0) {
        std::vector<int> nnz_parts;
        std::vector<int> m_parts;
        std::vector<int> groups;
        for(int k = 0; k < nbparts; k++) {
            nnz_parts.push_back(partitions[k].NumNonzeros());
            m_parts.push_back(partitions[k].dim(0));
        }

        //abcd::partitionWeights(groups, nnz_parts, parallel_cg);
        abcd::partitionWeights(groups, m_parts, parallel_cg);
        clog << "Groups : [" << groups[0] + 1;
        for(int k = 1; k < parallel_cg; k++) {
            clog  << ", "<< groups[k] - groups[k-1];
        }
        clog << "]" << endl;

        // first start by the master :
        for(int i = 0; i <= groups[0] ; i++)
            parts_id.push_back(i);

        for(int i = 1; i < parallel_cg ; i++) {
            // send to each CG-Master its starting and ending partitions ids
            std::vector<int> se;
            se.push_back(groups[i - 1]);
            se.push_back(groups[i]);

            inter_comm.send(i, 0, se);

            // send to each CG-Master the corresponding col_index and partitions
            for(int j = se[0] + 1; j <= se[1]; j++) {
                std::vector<int> sh;
                sh.push_back(partitions[j].dim(0));
                sh.push_back(partitions[j].dim(1));
                inter_comm.send(i, 1, partitions[j].NumNonzeros());
                inter_comm.send(i, 2, sh);
                inter_comm.send(i, 21, n);
                inter_comm.send(i, 3, partitions[j].colind_ptr(), partitions[j].NumNonzeros());
                inter_comm.send(i, 4, partitions[j].rowptr_ptr(), sh[0] + 1);
                inter_comm.send(i, 5, partitions[j].val_ptr(), partitions[j].NumNonzeros());

                inter_comm.send(i, 6, column_index[j]);

                if(icntl[10] > 0) inter_comm.send(i, 61, stC[j]);
            }
        }
        cout << "sent partitions" << endl;

        // Find the interconnections between the different CG_masters
        std::vector<std::vector<int> > group_column_index;
        int st = 0;
        for(int i = 0; i < parallel_cg ; i++) {
            std::vector<int> merge_index;
            merge_index.reserve(n);
            for(int j = st; j <= groups[i] ; j++) {
                std::copy(column_index[j].begin(), column_index[j].end(), back_inserter(merge_index));
            }
            std::sort(merge_index.begin(), merge_index.end());
            std::vector<int>::iterator last = std::unique(merge_index.begin(), merge_index.end());
            merge_index.erase(last, merge_index.end());
            group_column_index.push_back(merge_index);
            st = groups[i] + 1;

        }

        cout << "Merge done" << endl;

        // Send those interconnections to the other masters
        std::map<int, std::vector<int> > inter;
        for(int i = 0; i < parallel_cg; i++) {
            for(int j = i + 1; j < parallel_cg; j++) {
                std::vector<int> inter1;
                std::vector<int> inter2;
                std::vector<int>::iterator it1 = group_column_index[i].begin();
                std::vector<int>::iterator it2 = group_column_index[j].begin();

                inter1.reserve(group_column_index[j].size());
                inter2.reserve(group_column_index[j].size());

                while(it1 != group_column_index[i].end() && it2 != group_column_index[j].end()) {
                    if(*it1 < *it2) {
                        ++it1;
                    } else {
                        if(!(*it2 < *it1)) {
                            inter1.push_back(
				    it1 - group_column_index[i].begin()
                            );
                            inter2.push_back(
				    it2 - group_column_index[j].begin()
                            );

                            inter[i].push_back(
				    it1 - group_column_index[i].begin()
                            );
                            inter[j].push_back(
				    it1 - group_column_index[i].begin()
                            );
                        }
                        ++it2;
                    }
                }
                if(i == 0) {
                    col_interconnections[j] = inter1;
                } else {
                    inter_comm.send(i, 7, j);
                    inter_comm.send(i, 8, inter1);
                }
                inter_comm.send(j, 7, i);
                inter_comm.send(j, 8, inter2);
            }
        }


        for(int i = 1; i < parallel_cg; i++) {
            inter_comm.send(i, 7, -1);
            std::map<int, int> l_glob_to_local;
            for(int j = 0; j < group_column_index[i].size(); j++) {
                l_glob_to_local[group_column_index[i][j]] = j;
            }
        }

        //std::map<int, int> l_glob_to_local;
        //for(int j = 0; j < group_column_index[0].size(); j++) {
            //l_glob_to_local[group_column_index[0][j]] = j;
        //}

        //glob_to_local = l_glob_to_local;

        cout << "sent interconnections to others" << endl;


        if(groups[0] + 1 < nbparts){
            // Now that everybody got its partition, delete them from the master
            partitions.erase(partitions.begin() + groups[0] + 1, partitions.end());
            column_index.erase(column_index.begin() + groups[0] + 1, column_index.end());
            stC.erase(stC.begin() + groups[0] + 1, stC.end());
        }

        m_l = m;
        n_l = n;
        m = std::accumulate(partitions.begin(), partitions.end(), 0, sum_rows);
        nz = std::accumulate(partitions.begin(), partitions.end(), 0, sum_nnz);
    } else {
        std::vector<int> se;
        inter_comm.recv(0, 0, se);
        int sm = 0, snz = 0;
        for(int i = se[0] + 1; i <= se[1]; i++) {

            // The partition number to be added
            parts_id.push_back(i);

            // THe partition data
            std::vector<int> sh;
            int l_nz, l_m, l_n;
            inter_comm.recv(0, 1, l_nz);
            inter_comm.recv(0, 2, sh);
            inter_comm.recv(0, 21, n);
            l_m = sh[0];
            l_n = sh[1];

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
            if(icntl[10] > 0){
                inter_comm.recv(0, 61, stc);
                stC.push_back(stc);
            }

            // Create the matrix and push it in!
            CompRow_Mat_double mat(l_m, l_n, l_nz, l_v, l_irst, l_jcn);
            partitions.push_back(mat);

            sm += l_m;
            snz += l_nz;
        }

        while(true) {
            int the_other;
            std::vector<int> inter;
            inter_comm.recv(0, 7, the_other);
            if(the_other == -1) break;
            inter_comm.recv(0, 8, inter);
            col_interconnections[the_other] = inter;
        }


        // Set the number of rows and nnz handled by this CG Instance
        m = sm;
        nz = snz;
    }
    mpi::broadcast(inter_comm, m_l, 0);
    mpi::broadcast(inter_comm, n_l, 0);

    // Create a merge of the column indices
    std::vector<int> merge_index;
    for(int j = 0; j < partitions.size(); j++) {
        std::copy(column_index[j].begin(), column_index[j].end(), back_inserter(merge_index));
    }
    std::sort(merge_index.begin(), merge_index.end());
    std::vector<int>::iterator last = std::unique(merge_index.begin(), merge_index.end());
    merge_index.erase(last, merge_index.end());

    for(int j = 0; j < merge_index.size(); j++) {
        glob_to_local[merge_index[j]] = j;
    }

    for(int k = 0; k < partitions.size(); k++) {
        std::map<int, int> gt;
        std::map<int, int> ptg;
        for(int j = 0; j < column_index[k].size(); j++) {
            gt[column_index[k][j]] = j;
            ptg[j] = column_index[k][j];
        }
        part_to_glob.push_back(ptg);
        glob_to_part.push_back(gt);
    }

    // for each partition find a local column index for the previous merge
    std::vector<int> indices(partitions.size(), 0);

    local_column_index = std::vector<std::vector<int> >(partitions.size());

    for(int i = 0; i < partitions.size(); i++) {
        for(int j = 0; j < merge_index.size(); j++) {
            if(indices[i] >= column_index[i].size()) continue;
            if(column_index[i][indices[i]] == merge_index[j] &&
                    indices[i] < column_index[i].size() ) {

                local_column_index[i].push_back(j);
                indices[i]++;
            }
        }
    }

    n = merge_index.size();

    // Create a Communication map for ddot
    comm_map.assign(n, 1);
    for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
            it != col_interconnections.end(); it++) {
        // if I share data with it->first and I'm after him, let him compute!
        if(inter_comm.rank() > it->first) {
            for(std::vector<int>::iterator i = it->second.begin(); i != it->second.end(); i++) {
                if(comm_map[*i] == 1) comm_map[*i] = -1;
            }
        }
    }
    /*
     * ||A||_F
    nrmA = 0;
    for(int i=0; i<partitions.size(); i++)
        nrmA += squaredNorm(partitions[i]);

    mpi::all_reduce(inter_comm, &nrmA, 1,  &nrmMtx, std::plus<double>());
    nrmA = sqrt(nrmA);
    nrmMtx = sqrt(nrmMtx);
    */

    /* ||A||_inf */
    nrmA = 0;
    for(int i=0; i<partitions.size(); i++){
        for(int j=0; j<partitions[i].NumNonzeros(); j++){
            double cur = abs(partitions[i].val(j));
            if(cur > nrmA){
                nrmA = cur;
            }
        }
    }

    mpi::all_reduce(inter_comm, &nrmA, 1,  &nrmMtx, mpi::maximum<double>());
    //nrmA = sqrt(nrmA);
    //nrmMtx = sqrt(nrmMtx);
}

void abcd::distributeRhs()
{
    mpi::communicator world;

    mpi::broadcast(inter_comm, use_xf, 0);

    if(block_size < nrhs) block_size = nrhs;
    mpi::broadcast(inter_comm, block_size, 0);

    if(world.rank() == 0) {

        int r_pos = 0;
        // Build my part of the RHS
        int r = std::accumulate(partitions.begin(), partitions.end(), 0, sum_rows);

        if(rhs==NULL){
            //rhs = new double[r * nrhs];
            //for(int i=0; i<r*nrhs; i++) rhs[i] = ((double)rand()/(double)RAND_MAX);
            rhs = new double[n_l * nrhs];

            //for(int j = 0; j < obj.nrhs; j++){
                //for(int i = 0; i < obj.n_l; i++){
                    ////obj.rhs[i + j * obj.n_l] = j+1;
                    //obj.rhs[i + j * obj.n_l] = ((rand()%10)+j+1)/10; 
                //}
            //}
        }

        if(use_xf){

            nrmXf = 0;

            Xf = MV_ColMat_double(n_l, nrhs);

            for(int j = 0; j < nrhs; j++){
                VECTOR_double xf_col(n_l);
                for(int i = 0; i < n_l; i++) {
                    xf_col[i] = rhs[i + j * n_l];
                    if(abs(xf_col[i]) > nrmXf) nrmXf = abs(xf_col[i]);
                }
                Xf.setCol(xf_col, j);
            }

            B = smv(A, Xf);
            //

            // this is the augmented version!
            // ??? B should not change !
            //if(icntl[10]!=0){
                //for(int j = 0; j < B.size(); j++){
                    //VECTOR_double t(r, 0);
                    //t(MV_VecIndex(0, r-1)) = B[j](MV_VecIndex());
                    //B[j] = t;
                //}
            //}

            // for each rhs j
            for(int j = 0; j < B.dim(1); j++){
                for(int i = 0; i < B.dim(0); i++) {
                    rhs[i + j * B.dim(0)] = B(i,j);
                }
            }
        } else {
            B = MV_ColMat_double(m_l, block_size);

            Xf = MV_ColMat_double(A.dim(1), nrhs);
            for(int j = 0; j < nrhs; j++){
                for(int i = 0; i < A.dim(1); i++){
                    rhs[i + j * A.dim(1)] = j+1;
                    //rhs[i + j * n_l] = ((rand()%10)+j+1)/10; 
                }
            }

            nrmXf = 0;

            for(int j = 0; j < nrhs; j++){
                VECTOR_double xf_col(A.dim(1));
                for(int i = 0; i < A.dim(1); i++) {
                    xf_col[i] = rhs[i + j * A.dim(1)];
                    if(abs(xf_col[i]) > nrmXf) nrmXf = abs(xf_col[i]);
                }
                Xf.setCol(xf_col, j);
            }

            MV_ColMat_double BB = smv(A, Xf);

            for(int j = 0; j < nrhs; j++){
                //VECTOR_double t(rhs+j*m_l, m_l);
                //B.setCol(t, j);
                //B.push_back(t);
                B(MV_VecIndex(0, m_l-1), MV_VecIndex(0,nrhs-1)) = BB;
            }

            if(block_size > nrhs) {
                double *rdata = new double[m_l * (block_size - nrhs)];

                srand((unsigned)time(0)); 
                for(int i=0; i< m_l*(block_size-nrhs); i++){ 
                    rdata[i] = ((rand()%10)+1)/10; 
                    //rdata[i] = 2;
                }

                MV_ColMat_double RR(rdata, m_l, block_size-nrhs);
                B(MV_VecIndex(0,B.dim(0)-1),MV_VecIndex(nrhs,block_size-1)) = 
                    RR(MV_VecIndex(0,B.dim(0)-1), MV_VecIndex(0, block_size-nrhs - 1));
                delete[] rdata;
            }

        }

        r_pos += r;

        double *b_ptr = B.ptr();

        // for other masters except me!
        for(int k = 1; k < parallel_cg; k++) {
            // get the partitions that will be sent to the master
            int rows_for_k;
            inter_comm.recv(k, 16, rows_for_k);
            inter_comm.send(k, 17, nrhs);
            for(int j = 0; j < block_size; j++) {
                inter_comm.send(k, 18, b_ptr + r_pos + j * m_l, rows_for_k);
            }

            r_pos += rows_for_k;
        }

        if(!use_xf){
            MV_ColMat_double BB(m, block_size, 0);
            BB = B(MV_VecIndex(0, m-1), MV_VecIndex(0,block_size-1));
            B = MV_ColMat_double(m, block_size, 0);
            B = BB;
        }
    } else {
        inter_comm.send(0, 16, m);
        inter_comm.recv(0, 17, nrhs);

        //b = Eigen::MatrixXd(m, nrhs);
        //b.setZero();

        int size_rhs = m*block_size;
        rhs = new double[size_rhs];
        for(int i = 0; i < m * block_size; i++) rhs[i] = 0;

        B = MV_ColMat_double(m, block_size);
        for(int j = 0; j < block_size; j++) {

            inter_comm.recv(0, 18, rhs + j * m, m);
            VECTOR_double t(rhs+j*m, m);
            B.setCol(t, j);
        }
    }
    // and distribute max iterations
    mpi::broadcast(inter_comm, itmax, 0);
    mpi::broadcast(inter_comm, threshold, 0);
    mpi::broadcast(inter_comm, dcntl[10], 0);

    delete[] rhs;
}
