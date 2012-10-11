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
        std::vector<int> groups;
        for(int k = 0; k < nbparts; k++) {
            nnz_parts.push_back(parts[k].nonZeros());
        }

        abcd::partitionWeights(groups, nnz_parts, parallel_cg);

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
                sh.push_back(parts[j].rows());
                sh.push_back(parts[j].cols());
                inter_comm.send(i, 1, parts[j].nonZeros());
                inter_comm.send(i, 2, sh);
                inter_comm.send(i, 21, n);
                inter_comm.send(i, 3, parts[j].innerIndexPtr(), parts[j].nonZeros());
                inter_comm.send(i, 4, parts[j].outerIndexPtr(), sh[0] + 1);
                inter_comm.send(i, 5, parts[j].valuePtr(), parts[j].nonZeros());

                inter_comm.send(i, 6, column_index[j]);
            }
        }

        // Find the interconnections between the different CG_masters
        std::vector<std::vector<int> > group_column_index;
        int st = 0;
        for(int i = 0; i < parallel_cg ; i++) {
            std::vector<int> merge_index;
            for(int j = st; j <= groups[i] ; j++) {
                std::copy(column_index[j].begin(), column_index[j].end(), back_inserter(merge_index));
            }
            std::sort(merge_index.begin(), merge_index.end());
            std::vector<int>::iterator last = std::unique(merge_index.begin(), merge_index.end());
            merge_index.erase(last, merge_index.end());
            group_column_index.push_back(merge_index);
            st = groups[i] + 1;

        }

        std::map<int, std::vector<int> > inter;
        for(int i = 0; i < parallel_cg; i++) {
            for(int j = i + 1; j < parallel_cg; j++) {
                std::vector<int> inter1;
                std::vector<int> inter2;
                std::vector<int>::iterator it1 = group_column_index[i].begin();
                std::vector<int>::iterator it2 = group_column_index[j].begin();

                while(it1 != group_column_index[i].end() && it2 != group_column_index[j].end()) {
                    if(*it1 < *it2) {
                        ++it1;
                    } else {
                        if(!(*it2 < *it1)) {
                            inter1.push_back(it1 - group_column_index[i].begin());
                            inter2.push_back(it2 - group_column_index[j].begin());

                            inter[i].push_back(it1 - group_column_index[i].begin());
                            inter[j].push_back(it1 - group_column_index[i].begin());
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

        std::map<int, int> l_glob_to_local;
        for(int j = 0; j < group_column_index[0].size(); j++) {
            l_glob_to_local[group_column_index[0][j]] = j;
        }

        glob_to_local = l_glob_to_local;


        // Now that everybody got its partition, delete them from the master
        parts.erase(parts.begin() + groups[0] + 1, parts.end());
        column_index.erase(column_index.begin() + groups[0] + 1, column_index.end());
        local_column_index = column_index;

        m_l = m;
        m = std::accumulate(parts.begin(), parts.end(), 0, sum_rows);
        nz = std::accumulate(parts.begin(), parts.end(), 0, sum_nnz);
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
            int dec = ci[0];
            for(int i = 0; i < ci.size(); i++) ci[i] -= dec;
            local_column_index.push_back(ci);


            // Create the matrix and push it in!
            MappedSparseMatrix<double, RowMajor> mat(l_m, l_n, l_nz, l_irst, l_jcn, l_v);
            parts.push_back(SparseMatrix<double, RowMajor>(mat.middleRows(0, l_m)));

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

    std::vector<int> merge_index;
    for(int j = 0; j < parts.size(); j++) {
        std::copy(column_index[j].begin(), column_index[j].end(), back_inserter(merge_index));
    }
    std::sort(merge_index.begin(), merge_index.end());
    std::vector<int>::iterator last = std::unique(merge_index.begin(), merge_index.end());
    merge_index.erase(last, merge_index.end());

    int indices[parts.size()];
    for(int i = 0; i < parts.size(); i++) indices[i] = 0;
    local_column_index = std::vector<std::vector<int> >(parts.size());

    for(int j = 0; j < merge_index.size(); j++) {
        for(int i = 0; i < parts.size(); i++) {
            if(column_index[i][indices[i]] == merge_index[j]) {
                indices[i]++;
                local_column_index[i].push_back(j);
            }
        }
    }

    n = merge_index.size();

    // Create a Communication map for ddot
    comm_map = Eigen::VectorXi(n);
    comm_map.setOnes();
    for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
            it != col_interconnections.end(); it++) {
        if(inter_comm.rank() > it->first) {
            for(std::vector<int>::iterator i = it->second.begin(); i != it->second.end(); i++) {
                if(comm_map[*i] == 1) comm_map[*i] = -1;
            }
        }
    }
    nrmA = 0;
    for(int i=0; i<parts.size(); i++)
        nrmA += parts[i].squaredNorm();
    mpi::all_reduce(inter_comm, &nrmA, 1,  &nrmMtx, std::plus<double>());
    nrmA = sqrt(nrmA);
    nrmMtx = sqrt(nrmMtx);

}

void abcd::distributeRhs()
{
    mpi::communicator world;
    mpi::broadcast(inter_comm, block_size, 0);
    int s = std::max<int>(block_size, nrhs);
    if(s < 1) throw - 41;

    if(world.rank() == 0) {
        int r_pos = 0;
        // Build my part of the RHS
        int r = std::accumulate(parts.begin(), parts.end(), 0, sum_rows);
        b = Eigen::MatrixXd(r, s);

        if(s > nrhs) {
            Eigen::MatrixXd RR =  MatrixXd::Random(r, s - nrhs);
            //RR.setOnes();
            b.block(0, nrhs, r, s - nrhs) = RR;
        }

        for(int j = 0; j < nrhs; j++)
            for(int i = 0; i < r; i++) {
                b(i, j) = rhs[i + j * m];
            }

        r_pos += r;

        // for other masters except me!
        for(int k = 1; k < parallel_cg; k++) {
            // get the partitions that will be sent to the master
            int rows_for_k;
            inter_comm.recv(k, 16, rows_for_k);
            inter_comm.send(k, 17, nrhs);
            for(int j = 0; j < nrhs; j++) {
                inter_comm.send(k, 18, rhs + r_pos + j * m, rows_for_k);
            }

            r_pos += rows_for_k;
        }
    } else {
        inter_comm.send(0, 16, m);
        inter_comm.recv(0, 17, nrhs);

        b = Eigen::MatrixXd(m, s);
        b.setZero();

        rhs = new double[m * nrhs];
        for(int i = 0; i < m * nrhs; i++) rhs[i] = 0;

        for(int j = 0; j < nrhs; j++) {

            inter_comm.recv(0, 18, rhs + j * m, m);

            for(int i = 0; i < m; i++) {
                b(i, j) = rhs[i + j * m];
            }

        }

        if(s > nrhs) {
            Eigen::MatrixXd RR =  MatrixXd::Random(m, s - nrhs);
            //RR.setOnes();
            b.block(0, nrhs, m, s - nrhs) = RR;
        }
    }
    
    // and distribute max iterations
    mpi::broadcast(inter_comm, itmax, 0);
}






















