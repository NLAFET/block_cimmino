#include "abcd.h"

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
            std::vector<int> se = {groups[i - 1], groups[i]};
            inter_comm.send(i, 0, se);

            // send to each CG-Master the corresponding col_index and partitions
            for(int j = se[0] + 1; j <= se[1]; j++) {
                std::vector<int> sh;
                sh.push_back(parts[j].rows());
                sh.push_back(parts[j].cols());
                inter_comm.send(i, 1, parts[j].nonZeros());
                inter_comm.send(i, 2, sh);
                inter_comm.send(i, 3, parts[j].innerIndexPtr(), parts[j].nonZeros());
                inter_comm.send(i, 4, parts[j].outerIndexPtr(), sh[0] + 1);
                inter_comm.send(i, 5, parts[j].valuePtr(), parts[j].nonZeros());

                inter_comm.send(i, 6, columns_index[j]);
            }
        }

        // Find the interconnections between the different CG_masters
        std::vector<std::vector<int> > group_column_index;
        int st = 0;
        for(int i = 0; i < parallel_cg ; i++) {
            for(int j = st; j <= groups[i] ; j++) {
                std::vector<int> merge_index;
                std::copy(columns_index[j].begin(), columns_index[j].end(), back_inserter(merge_index));
                std::sort(merge_index.begin(), merge_index.end());
                std::vector<int>::iterator last = std::unique(merge_index.begin(), merge_index.end());
                merge_index.erase(last, merge_index.end());
                group_column_index.push_back(merge_index);
            }
            st = groups[i] + 1;
        }
        for(int i = 0; i < parallel_cg; i++) {
            for(int j = i + 1; j < parallel_cg; j++) {
                std::vector<int> inter;
                set_intersection(group_column_index[i].begin(), group_column_index[i].end(),
                                 group_column_index[j].begin(), group_column_index[j].end(),
                                 back_inserter(inter));
                if(i == 0) {
                    col_interconnections[j] = inter;
                } else {
                    inter_comm.send(i, 7, j);
                    inter_comm.send(i, 8, inter);
                }
                inter_comm.send(j, 7, i);
                inter_comm.send(j, 8, inter);
            }
        }


        // Now that everybody got its partition, delete them from the master
        parts.erase(parts.begin() + groups[0] + 1, parts.end());
        columns_index.erase(columns_index.begin() + groups[0] + 1, columns_index.end());
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
            columns_index.push_back(ci);


            // Create the matrix and push it in!
            MappedSparseMatrix<double, RowMajor> mat(l_m, l_n, l_nz, l_irst, l_jcn, l_v);
            parts.push_back(SparseMatrix<double, RowMajor>(mat.middleRows(0, l_m)));

            sm += l_m;
            snz += l_nz;
        }

        for(int i = inter_comm.rank(); i < parallel_cg; i++) {
            int the_other;
            std::vector<int> inter;
            inter_comm.recv(0, 7, the_other);
            inter_comm.recv(0, 8, inter);
            col_interconnections[the_other] = inter;
        }

        // Set the number of rows and nnz handled by this CG Instance
        m = sm;
        nz = snz;
    }
}


