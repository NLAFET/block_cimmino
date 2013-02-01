#include <abcd.h>
#include <Eigen/src/Geometry/Translation.h>
void abcd::initializeMumps()
{
    initializeMumps(false);
}

void abcd::initializeMumps(bool local)
{
    mpi::communicator world;
    // The first run of MUMPS is local to CG-masters
    std::vector<int> r;
    if(local) {
        r.push_back(inter_comm.rank());
        mpi::group grp = inter_comm.group().include(r.begin(), r.end());
        intra_comm = mpi::communicator(inter_comm, grp);
    } else {
        if(instance_type == 0) {
            r.push_back(world.rank());
        } else {
            r.push_back(my_master);
        }
        std::copy(my_slaves.begin(), my_slaves.end(), std::back_inserter(r));

        mpi::group grp = world.group().include(r.begin(), r.end());
        intra_comm = mpi::communicator(world, grp);
    }

    mumps.sym = 2;
    mumps.par = 1;
    mumps.job = -1;
    mumps.comm_fortran = MPI_Comm_c2f((MPI_Comm) intra_comm);

    dmumps_c(&mumps);
    if(getMumpsInfo(1) != 0) throw - 100 + getMumpsInfo(1);

    setMumpsIcntl(1, -1);
    setMumpsIcntl(2, -1);
    setMumpsIcntl(3, -1);

    setMumpsIcntl(14, 50);
    setMumpsIcntl(12, 2);
    setMumpsIcntl(6, 5);
    setMumpsIcntl(7, 5);
    setMumpsIcntl(8, -2);
    setMumpsIcntl(27, 16);
}

void abcd::createAugmentedSystems()
{
    SparseMatrix<double, RowMajor> G;
    
    m_n = 0;
    m_nz = 0;
    // for performance, compute total nnz and size of the matrix
    for(int j = 0; j < partitions.size(); j++) {
        m_n += partitions[j].dim(0) + partitions[j].dim(1);
        m_nz += partitions[j].dim(1) + partitions[j].NumNonzeros();
    }

    // Allocate the data for mumps
    mumps.n = m_n;
    mumps.nz = m_nz;
    mumps.irn = new int[m_nz];
    mumps.jcn = new int[m_nz];
    mumps.a = new double[m_nz];

    // Use Fortran array => start from 1
    int i_pos = 1;
    int j_pos = 1;
    int st = 0;

    for(int p = 0; p < partitions.size(); p++) {

        // fill the identity
        for(int i = 0; i < partitions[p].dim(1); i++) {
            mumps.irn[st + i] = i_pos + i;
            mumps.jcn[st + i] = j_pos + i;
            mumps.a[st + i] = 1;
        }

        // we get down by nb_cols
        i_pos += partitions[p].dim(1);
        // we added nb_cols elements
        st += partitions[p].dim(1);

        for(int k = 0; k < partitions[p].dim(0); k++) {
            for(int j = partitions[p].row_ptr(k); j < partitions[p].row_ptr(k + 1); j++) {
                mumps.irn[st] = i_pos + k;
                mumps.jcn[st] = j_pos + partitions[p].col_ind(j);
                mumps.a[st] = partitions[p].val(j);
                st++;
            }
        }

        i_pos += partitions[p].dim(0);
        j_pos += partitions[p].dim(1) + partitions[p].dim(0);

    }
//     The data given to MUMPS
//     for(int k = 0; k < st; k++) {
//         cout << mumps.irn[k] << " ";
//         cout << mumps.jcn[k] << " ";
//         cout << mumps.a[k] << endl;
//     }
//
//     cout << l_nnz << " " << st << endl;
//     cout << l_n<< " " << i_pos << endl;
}

void abcd::allocateMumpsSlaves()
{
    mpi::communicator world;

    if(instance_type == 0) {
        std::vector<int> flops(inter_comm.size());
        std::vector<dipair> flops_s;
        std::vector<int> slaves_for_me(inter_comm.size());
        std::vector<int> slaves_for_me_t(inter_comm.size());
        mpi::all_gather(inter_comm, (int) getMumpsRinfo(1), flops);

        for(int idx = 0; idx < inter_comm.size(); idx++) {
            flops_s.push_back(dipair(flops[idx], idx));
        }

        std::sort(flops_s.begin(), flops_s.end(), ip_comp);


        int s = std::accumulate(flops.begin(), flops.end(), 0);
        std::vector<double> shares;
        int nb_slaves = world.size() - inter_comm.size();
        int slaves_left = nb_slaves;
        double top = 1, low = 0.90;

        for(int i = 0; i < inter_comm.size() && slaves_left > 0 ; i++) {
            shares.push_back(
		    ((double)flops_s[i].first / (double) s) * nb_slaves
            );
        }

        while(slaves_left > 0) {
            for(int i = 0; i < inter_comm.size() && slaves_left > 0 ; i++) {
                int share_of_slaves = 0;
                if(shares[i] < 0) continue;
                if((shares[i] - floor(shares[i])) >= low  &&
                        (shares[i] - floor(shares[i])) < top) {
		    share_of_slaves = ceil(shares[i]) < slaves_left ?
                        ceil(shares[i]) : slaves_left;
                } else {
                    share_of_slaves = floor(shares[i]) < slaves_left ?
                        floor(shares[i]) : slaves_left;
                }
                slaves_for_me_t[i] += share_of_slaves;
                slaves_left -= share_of_slaves;
                shares[i] -= share_of_slaves;
            }
            top -= 0.10;
            low -= 0.10;
        }
        for(int idx = 0; idx < inter_comm.size(); idx++) {
            slaves_for_me[ flops_s[idx].second ] = slaves_for_me_t[idx];
        }

        int current_slave = inter_comm.size();
        for(int your_master = 0 ; your_master < inter_comm.size();
                your_master++) {

            for(int i = 0; i < slaves_for_me[your_master]; i++) {

                // let the master handle this!
                if(inter_comm.rank() == 0)
                    world.send(current_slave, 11,  your_master);

                if(inter_comm.rank() == your_master) {
                    my_slaves.push_back(current_slave);
                }

                current_slave++;
            }
        }
        // Now that the slaves know who's their daddy,
        //tell who are their brothers
        for(std::vector<int>::iterator slave = my_slaves.begin();
                slave != my_slaves.end(); slave++) {
            world.send(*slave, 12, my_slaves);
        }

    } else {
        world.recv(0, 11, my_master);
        //Store in my_slaves my brothers in slavery
        world.recv(my_master, 12, my_slaves);
    }
}

void abcd::analyseAugmentedSystems()
{
    mumps.job = 1;

    double t = MPI_Wtime();
    //cout << "mumps.n " << mumps.n << endl;
    //cout << "mumps.nz " << mumps.nz << endl;
    dmumps_c(&mumps);
    t = MPI_Wtime() - t;

    if(getMumpsInfo(1) != 0) throw - 10000 + 100 * getMumpsInfo(1) - getMumpsInfo(2);

    if(instance_type == 0) {
        double flop = getMumpsRinfo(1);
        int prec = cout.precision();
        cout.precision(2);
        cout << string(32, '-') << endl
             << "| MUMPS ANALYSIS on MA " << setw(7) << inter_comm.rank() << " |" << endl
             << string(32, '-') << endl
             << "| Flops estimate: " << setw(6) << scientific << flop << string(4, ' ') << " |" << endl
             << "| Time:           " << setw(6) << t << " sec |" << endl
             << string(32, '-') << endl;
        cout.precision(prec);
    }
}


void abcd::factorizeAugmentedSystems()
{
    mpi::communicator world;
    mumps.job = 2;

    double t = MPI_Wtime();
    dmumps_c(&mumps);
    t = MPI_Wtime() - t;

    if(getMumpsInfo(1) != 0) throw - 10000 + 100 * getMumpsInfo(1) - getMumpsInfo(2);

    if(instance_type == 0) {
        double flop = getMumpsRinfo(1);
        int prec = cout.precision();
        cout.precision(2);
        cout << string(32, '-') << endl
             << "| MUMPS FACTORIZ on MA " << setw(7) << inter_comm.rank() << " |" << endl
             << string(32, '-') << endl
             << "| Flops estimate: " << setw(6) << scientific << flop << string(4, ' ') << " |" << endl
             << "| Time:           " << setw(6) << t << " sec |" << endl
             << string(32, '-') << endl;
        cout.precision(prec);
    }
}

MV_ColMat_double abcd::sumProject(double alpha, MV_ColMat_double &Rhs, double beta, MV_ColMat_double &X)
{
    mpi::communicator world;
    int s = X.dim(1);
    // Build the mumps rhs
    mumps.rhs = new double[mumps.n * s];
    int pos = 0;
    int b_pos = 0;
    for(int k = 0; k < partitions.size(); k++) {
        MV_ColMat_double r(partitions[k].dim(0), s, 0);
        MV_ColMat_double compressed_x(partitions[k].dim(1), s, 0);

        int x_pos = 0;
        for(int i = 0; i < local_column_index[k].size(); i++) {
            int ci = local_column_index[k][i];
            for(int j = 0; j < s; j++) {
                assert(x_pos < compressed_x.dim(0));
                assert(ci < X.dim(0));
                compressed_x(x_pos, j) = X(ci, j);
            }
            x_pos++;
        }

        // avoid useless operations
        if(beta != 0){
            r = smv(partitions[k], compressed_x) * beta;
        }
        if(alpha != 0){
            //r = r + alpha * B.block(b_pos, 0, parts[k].rows(), s);
            r = r + B(MV_VecIndex(b_pos, b_pos + partitions[k].dim(0) - 1),
                              MV_VecIndex(0, s -1)) * alpha;
        }
        
        b_pos += partitions[k].dim(0);

        for(int r_p = 0; r_p < s; r_p++) {
            for(int i = pos; i < pos + partitions[k].dim(1); i++) {
                mumps.rhs[i + r_p * mumps.n] = 0;
            }
            int j = 0;
            for(int i = pos + partitions[k].dim(1);
                    i < pos + partitions[k].dim(1) + partitions[k].dim(0); i++) {
                mumps.rhs[i + r_p * mumps.n] = r(j++, r_p);
            }
        }

        pos += partitions[k].dim(1) + partitions[k].dim(0);

    }
    mumps.nrhs = s;
    mumps.lrhs = mumps.n;
    mumps.job = 3;
    double t = MPI_Wtime();
    dmumps_c(&mumps);
    t = MPI_Wtime() - t;
    cout << "[" << inter_comm.rank() << "] Time spent in direct solver : "
        << t << endl;

    MV_ColMat_double Sol(mumps.rhs, mumps.n, s);
    MV_ColMat_double Delta(X.dim(0), X.dim(1), 0);

    int x_pos = 0;
    for(int k = 0; k < partitions.size(); k++) {
        for(int i = 0; i < local_column_index[k].size(); i++) {
            int ci = local_column_index[k][i];
            for(int j = 0; j < s; j++) {
                Delta(ci, j) = Delta(ci, j) + Sol(x_pos, j) ;
            }
            x_pos++;
        }
        x_pos += partitions[k].dim(0);
    }
    // Where the other Deltas are going to be summed
    MV_ColMat_double Others(X.dim(0), X.dim(1), 0);


    for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin(); it != col_interconnections.end(); it++) {

        // Prepare the data to be sent
        std::vector<double> itc, otc;
        for(int j = 0; j < s; j++) {
            for(std::vector<int>::iterator i = it->second.begin(); i != it->second.end(); i++) {
                itc.push_back(Delta(*i, j));
            }
        }

        // EVEN -> ODD && ODD -> EVEN
        if(inter_comm.rank() % 2 == 0) {
            if(it->first % 2 != 0) {
                inter_comm.send(it->first, 31, itc);
                inter_comm.recv(it->first, 32, otc);
            } else { // EVEN -> EVEN
                if(it->first > inter_comm.rank()) {
                    inter_comm.send(it->first, 33, itc);
                    inter_comm.recv(it->first, 34, otc);
                } else {
                    int z;
                    inter_comm.recv(it->first, 33, otc);
                    inter_comm.send(it->first, 34, itc);
                }
            }
        } else {
            if(it->first % 2 == 0) {
                inter_comm.recv(it->first, 31, otc);
                inter_comm.send(it->first, 32, itc);
            } else { // ODD -> ODD
                if(it->first > inter_comm.rank()) {
                    inter_comm.send(it->first, 33, itc);
                    inter_comm.recv(it->first, 34, otc);
                } else {
                    inter_comm.recv(it->first, 33, otc);
                    inter_comm.send(it->first, 34, itc);
                }
            }
        }

        // Uncompress data and sum it inside Others
        int p = 0;
        for(int j = 0; j < s; j++) {
            for(std::vector<int>::iterator i = it->second.begin(); i != it->second.end(); i++) {
                Others(*i, j) += otc[p++];
            }
        }

    }

    // Now sum the data to Delta
    Delta += Others;

    return Delta;
}


