#include <abcd.h>
#include <vect_utils.h>

void abcd::distributeData()
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

        abcd::partitioning(p_sets, m_parts, parallel_cg);

        clog << "Groups : " ;
        for(int k = 0; k < parallel_cg; k++) {
            clog << "{";
            for(int j = 0; j < p_sets[k].size() - 1; j++)
                clog << p_sets[k][j] << ", ";
            clog << p_sets[k][p_sets[k].size() -1 ];
            clog << "} ";
            //clog  << ", "<< p_sets[k].size();
        }
        clog << endl;

        for(int i = 1; i < parallel_cg ; i++) {
            inter_comm.send(i, 0, p_sets[i]);

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

        m_l = m;
        n_l = n;
        m = 0;
        nz = 0;

        if(parallel_cg != 1){
            std::vector<std::vector<int> > cis;
            std::vector<int> stcs;

            for(int i = 0; i < p_sets[0].size(); i++){
                int j = p_sets[0][i];
                partitions.push_back(parts[j]);
                m += parts[j].dim(0);
                nz += parts[j].NumNonzeros();

                cis.push_back(column_index[j]);
                if(icntl[10] != 0) stcs.push_back(stC[j]);
            }
            parts.clear();
            column_index.clear();
            nb_local_parts = partitions.size();
            if(icntl[10] != 0) stC.clear();
            for(int i = 0; i < nb_local_parts; i++){
                //partitions[i] = tp[i];
                column_index.push_back(cis[i]);
                if(icntl[10] != 0) stC.push_back(stcs[i]);
            }
        } else {
            for(int i = 0; i < parts.size(); i++){
                partitions.push_back(parts[i]);
                m += parts[i].dim(0);
                nz += parts[i].NumNonzeros();
            }
            parts.clear();
            nb_local_parts = partitions.size();
        }


    } else {
        std::vector<int> se;
        inter_comm.recv(0, 0, se);
        int sm = 0, snz = 0;
        for(int i = 0; i < se.size(); i++) {


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
        nb_local_parts = partitions.size();


        // Set the number of rows and nnz handled by this CG Instance
        m = sm;
        nz = snz;
    }
    mpi::broadcast(inter_comm, m_l, 0);
    mpi::broadcast(inter_comm, n_l, 0);
}

void abcd::createInterconnections()
{
    // we need the merge of column indices 
    std::vector<int> merge_index = mergeSortedVectors(column_index);

    // for ABCD, we need a global to local indices so that we can
    // identify which column in C is linked to 
    for(int j = 0; j < merge_index.size(); j++) {
        glob_to_local[merge_index[j]] = j;
        glob_to_local_ind.push_back(merge_index[j]);
    }
    // defines the starting point of C in the local columns
    st_c_part_it = glob_to_local_ind.end();


    while(*(st_c_part_it - 1) >= n_o) st_c_part_it--;

    st_c_part = st_c_part_it - glob_to_local_ind.begin();

    // for each partition find a local column index for the previous merge
    local_column_index = std::vector<std::vector<int> >(nb_local_parts);

    for(int p = 0; p < nb_local_parts; p++) {
        int i = 0, j = 0;
        while (i < column_index[p].size() && j < merge_index.size()) {

            if(column_index[p][i] == merge_index[j]) {
                local_column_index[p].push_back(j);
                i++;
            }

            j++;
        }
    }


    cout << "Merge done" << endl;

    // first build partitions to local vector
}
