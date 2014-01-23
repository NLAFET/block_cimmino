#include <abcd.h>
#include <vect_utils.h>

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
            nnz_parts.push_back(parts[k].NumNonzeros());
            m_parts.push_back(parts[k].dim(0));
        }

        //abcd::partitionWeights(groups, nnz_parts, parallel_cg);
        //abcd::partitionWeights(groups, m_parts, parallel_cg);
        abcd::partitioning(p_sets, m_parts, parallel_cg);
        //exit(0);

        clog << "Groups : [" << p_sets[0].size();
        for(int k = 1; k < parallel_cg; k++) {
            clog  << ", "<< p_sets[k].size();
        }
        clog << "]" << endl;

        // first start by the master :
        //for(int i = 0; i <= groups[0] ; i++)
            //parts_id.push_back(i);

        for(int i = 1; i < parallel_cg ; i++) {
            // send to each CG-Master its starting and ending partitions ids
            //std::vector<int> se;
            //se.push_back(groups[i - 1]);
            //se.push_back(groups[i]);

            //inter_comm.send(i, 0, se);
            inter_comm.send(i, 0, p_sets[i]);

            // send to each CG-Master the corresponding col_index and partitions
            //for(int j = se[0] + 1; j <= se[1]; j++) {
            for(int k = 0; k < p_sets[i].size(); k++){
                int j = p_sets[i][k];

                std::vector<int> sh;
                sh.push_back(parts[j].dim(0));
                sh.push_back(parts[j].dim(1));
                inter_comm.send(i, 1, parts[j].NumNonzeros());
                inter_comm.send(i, 2, sh);
                inter_comm.send(i, 21, n);
                inter_comm.send(i, 3, parts[j].colind_ptr(), parts[j].NumNonzeros());
                inter_comm.send(i, 4, parts[j].rowptr_ptr(), sh[0] + 1);
                inter_comm.send(i, 5, parts[j].val_ptr(), parts[j].NumNonzeros());

                inter_comm.send(i, 6, column_index[j]);


                if(icntl[10] > 0) inter_comm.send(i, 61, stC[j]);
            }

        }
        cout << "sent partitions" << endl;

        // Find the interconnections between the different CG_masters
        std::vector<std::vector<int> > group_column_index;
        int st = 0;
        for(int i = 0; i < parallel_cg ; i++) {
            std::vector<std::vector<int> > cis;
            for(int k = 0; k < p_sets[i].size(); k++){
                int j = p_sets[i][k];
                cis.push_back(column_index[j]);
            }
            std::vector<int> merge_index = mergeSortedVectors(cis);
            group_column_index.push_back(merge_index);
        }

        cout << "Merge done" << endl;

        mpi::broadcast(inter_comm, group_column_index, 0);

        // Send those interconnections to the other masters
        std::map<int, std::vector<int> > inter;
        //for(int i = 0; i < parallel_cg; i++) 
        {
            int i = 0;
            for(int j = i + 1; j < parallel_cg; j++) {
                std::pair<std::vector<int>, std::vector<int> > p =
                    getIntersectionIndices(
                            group_column_index[i], group_column_index[j]
                            );
                col_interconnections[j] = p.first;
            }
        }

        cout << "sent interconnections to others" << endl;

        if(parallel_cg != 1){
            //std::map<int, CompRow_Mat_double> tp;
            std::vector<std::vector<int> > cis;
            std::vector<int> stcs;

            for(int i = 0; i < p_sets[0].size(); i++){
                int j = p_sets[0][i];
                //tp[i] = partitions[j];
                partitions.push_back(parts[j]);
                cis.push_back(column_index[j]);
                if(icntl[10] != 0) stcs.push_back(stC[j]);
            }
            parts.clear();
            column_index.clear();
            if(icntl[10] != 0) stC.clear();
            for(int i = 0; i < partitions.size(); i++){
                //partitions[i] = tp[i];
                column_index.push_back(cis[i]);
                if(icntl[10] != 0) stC.push_back(stcs[i]);
            }
        } else {
            for(int i = 0; i < parts.size(); i++){
                partitions.push_back(parts[i]);
            }
            parts.clear();
        }

        m_l = m;
        n_l = n;
        m = 0;
        nz = 0;
        for(int i = 0; i < partitions.size(); i++){
            m += partitions[i].dim(0);
            nz += partitions[i].NumNonzeros();
        }
        //m = std::accumulate(partitions.begin(), partitions.end(), 0, sum_rows);
        //nz = std::accumulate(partitions.begin(), partitions.end(), 0, sum_nnz);
    } else {
        std::vector<int> se;
        inter_comm.recv(0, 0, se);
        int sm = 0, snz = 0;
        //for(int i = se[0] + 1; i <= se[1]; i++) {
        for(int i = 0; i < se.size(); i++) {

            // The partition number to be added
            //parts_id.push_back(i);

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
            delete [] l_jcn, l_irst, l_v;
        }

        // Find the interconnections between the different CG_masters
        std::vector<std::vector<int> > group_column_index;
        mpi::broadcast(inter_comm, group_column_index, 0);

        {
            int i = inter_comm.rank();
            for(int j = 0; j < group_column_index.size(); j++) {
                if(j == inter_comm.rank()) continue;
                std::pair<std::vector<int>, std::vector<int> > p =
                    getIntersectionIndices(
                            group_column_index[i], group_column_index[j]
                            );
                col_interconnections[j] = p.first;
            }
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
        glob_to_local_ind.push_back(merge_index[j]);
    }
    st_c_part_it = glob_to_local_ind.end();

    while(*(st_c_part_it - 1) >= n_o) st_c_part_it--;

    st_c_part = st_c_part_it - glob_to_local_ind.begin();

    for(std::vector<int>::iterator iti = st_c_part_it; iti != glob_to_local_ind.end(); iti++){
        glob_to_local_c[*iti - n_o] = iti - glob_to_local_ind.begin();
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

    // test!
    fast_local_column_index = new int *[partitions.size()];
    for(int i = 0; i < partitions.size(); i++) {
        indices[i] = 0;
        int p = 0;
        fast_local_column_index[i] = new int[local_column_index[i].size()];
        for(int j = 0; j < merge_index.size(); j++) {
            if(indices[i] >= column_index[i].size()) continue;
            if(column_index[i][indices[i]] == merge_index[j] &&
                    indices[i] < column_index[i].size() ) {

                fast_local_column_index[i][p] = j;
                p++;
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

    std::map<int, std::vector<int>::iterator > debut;
    std::map<int, int > their_job;
    {
        std::map<int, std::vector<int> > their_cols;
        std::vector<mpi::request> reqs_c;
        for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
                it != col_interconnections.end(); it++) {
            reqs_c.push_back(inter_comm.irecv(it->first, 40, their_cols[it->first]));
            reqs_c.push_back(inter_comm.isend(it->first, 40, glob_to_local_ind));
        }
        mpi::wait_all(reqs_c.begin(), reqs_c.end());
        reqs_c.clear();

        for(std::map<int, std::vector<int> >::iterator it = their_cols.begin();
                it != their_cols.end(); it++) {
            their_job[it->first] = it->second.size();
            set_intersection(
                    glob_to_local_ind.begin(), glob_to_local_ind.end(),
                    it->second.begin(), it->second.end(),
                    back_inserter(col_inter[it->first])
                    );
        }

    }

    mpi::broadcast(inter_comm, selected_S_columns, 0);
    mpi::broadcast(inter_comm, skipped_S_columns, 0);

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
        //int r = std::accumulate(partitions.begin(), partitions.end(), 0, sum_rows);
        int r = 0;
        for(int i = 0; i < partitions.size(); i++){
            r += partitions[i].dim(0);
        }

        if(rhs==NULL){
            rhs = new double[n_l * nrhs];

            srand(10); 
            B = MV_ColMat_double(m_l, block_size);

            nrmXf = 0;
            Xf = MV_ColMat_double(A.dim(1), nrhs);
            for(int j = 0; j < nrhs; j++){
                for(int i = 0; i < A.dim(1); i++){
                    //rhs[i + j * A.dim(1)] = 1 / dcol_[i];
                    //rhs[i + j * A.dim(1)] = 1;
                    //rhs[i + j * A.dim(1)] = j+1;
                    //rhs[i + j * n_l] = ((rand()%n_l)+j+1)/((double) n_l); 
                    rhs[i + j * n_l] = (double)((rand())%100+1)/99.0;

                    if(nrmXf < abs(rhs[i + j * n_l])) nrmXf = abs(rhs[i + j * n_l]);
                }
            }

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
        } else {
            B = MV_ColMat_double(m_l, block_size, 0);
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

        if(block_size > nrhs) {
            double *rdata = new double[m_l * (block_size - nrhs)];

            srand(n_l); 
            for(int i=0; i< m_l*(block_size-nrhs); i++){ 
                rdata[i] = (double)((rand())%10)/99.9 + 1;
                //rdata[i] = i+1;
            }
            MV_ColMat_double BR(rdata, m_l, block_size - nrhs, MV_Matrix_::ref);
            MV_ColMat_double RR = smv(A, BR);

            B(MV_VecIndex(0,B.dim(0)-1),MV_VecIndex(nrhs,block_size-1)) = 
                RR(MV_VecIndex(0,B.dim(0)-1), MV_VecIndex(0, block_size-nrhs - 1));
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
            for(int j = 0; j < block_size; j++) {
                for(int i = 0; i < p_sets[k].size(); i++){
                    int p = p_sets[k][i];
                    inter_comm.send(k, 18, b_ptr + strow[p] + j * m_l, nbrows[p]);
                }
            }
        }

        if(!use_xf){
            MV_ColMat_double BB(B);

            B = MV_ColMat_double(m, block_size, 0);
            int pos = 0;
            for(int i = 0; i < p_sets[0].size(); i++){
                int p = p_sets[0][i];
                for(int j = 0; j < nbrows[p]; j++){
                    for(int k = 0; k < block_size; k++){
                        B(pos, k) = BB(strow[p] + j, k);
                    }
                    pos++;
                }
            }
        }
    } else {
        Xf = MV_ColMat_double(n_o, 1, 0);
        double *xf_ptr = Xf.ptr();
        inter_comm.recv(0, 171, xf_ptr, n_o);

        inter_comm.recv(0, 17, nrhs);

        //b = Eigen::MatrixXd(m, nrhs);
        //b.setZero();

        int size_rhs = m*block_size;
        rhs = new double[size_rhs];
        for(int i = 0; i < m * block_size; i++) rhs[i] = 0;

        B = MV_ColMat_double(m, block_size, 0);
        for(int j = 0; j < block_size; j++) {
            int p = 0;
            for(int k = 0; k < partitions.size(); k++){
                inter_comm.recv(0, 18, rhs + p + j * m, partitions[k].dim(0));
                p+= partitions[k].dim(0);
            }

            VECTOR_double t(rhs+j*m, m);
            B.setCol(t, j);
        }
    }
    // and distribute max iterations
    mpi::broadcast(inter_comm, itmax, 0);
    mpi::broadcast(inter_comm, threshold, 0);
    mpi::broadcast(inter_comm, dcntl[10], 0);

    delete[] rhs;
    A = CompRow_Mat_double();
}

void abcd::distributeNewRhs()
{
    mpi::communicator world;


    if(world.rank() == 0) {

        int r_pos = 0;
        // Build my part of the RHS
        //int r = std::accumulate(partitions.begin(), partitions.end(), 0, sum_rows);
        int r = 0;
        for(int i = 0; i < partitions.size(); i++){
            r += partitions[i].dim(0);
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

    delete[] rhs;
}
