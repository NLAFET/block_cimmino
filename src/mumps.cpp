#include <abcd.h>

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

    mumps.sym = 0;
    mumps.par = 1;
    mumps.job = -1;
    mumps.comm_fortran = MPI_Comm_c2f((MPI_Comm) intra_comm);

    dmumps_c(&mumps);
    //cout << world.rank() << " << >> " << getMumpsInfo(1) << endl;
    if(getMumpsInfo(1) != 0) throw getMumpsInfo(1);

    setMumpsIcntl(1, -1);
    setMumpsIcntl(2, -1);
    setMumpsIcntl(3, -1);

    //if(inter_comm.rank() == 0){ 
        //strcpy(mumps.write_problem, "/scratch/mzenadi/pb6.mtx");
        //setMumpsIcntl(1, 6);
        //setMumpsIcntl(2, 6);
        //setMumpsIcntl(3, 6);
    //}


    setMumpsIcntl(6, 5);
    setMumpsIcntl(7, 5);
    setMumpsIcntl(8, -2);
    //setMumpsIcntl(11, 1);
    //setMumpsIcntl(12, 2);
    setMumpsIcntl(14, 50);
}

void abcd::createAugmentedSystems()
{
    m_n = 0;
    m_nz = 0;
    // for performance, compute total nnz and size of the matrix
    for(int j = 0; j < partitions.size(); j++) {
        m_n += partitions[j].dim(0) + partitions[j].dim(1);
        //m_nz += partitions[j].dim(1) + partitions[j].NumNonzeros();
        m_nz += partitions[j].dim(1) + 2 * partitions[j].NumNonzeros();
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
                mumps.jcn[st] = i_pos + k;
                mumps.irn[st] = j_pos + partitions[p].col_ind(j);
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
        std::vector<long> flops(inter_comm.size());
        std::vector<std::pair<long, int> > flops_s;
        std::vector<int> slaves_for_me(inter_comm.size());
        std::vector<int> slaves_for_me_t(inter_comm.size());
        mpi::all_gather(inter_comm, (long) getMumpsRinfo(1), flops);

        long s = 0;
        for(int idx = 0; idx < inter_comm.size(); idx++) {
            flops_s.push_back(std::pair<long,int>(flops[idx], idx));
            s+=flops[idx];
        }

        std::sort(flops_s.begin(), flops_s.end(), ip_comp);

        std::vector<double> shares;
        int nb_slaves = world.size() - inter_comm.size();
        int slaves_left = nb_slaves;
        double top = 1, low = 0.90;

        for(int i = 0; i < inter_comm.size() && slaves_left > 0 ; i++) {
            shares.push_back( ((double)flops_s[i].first / (double) s) * nb_slaves);
        }

        while(slaves_left > 0) {
            for(int i = 0; i < inter_comm.size() && slaves_left > 0 ; i++) {
                int share_of_slaves = 0;
                if(shares[i] < 0) continue;
                if((shares[i] - floor(shares[i])) >= low  &&
                        (shares[i] - floor(shares[i])) < top) {
                    share_of_slaves = ceil(shares[i]) < slaves_left ?  ceil(shares[i]) : slaves_left;
                } else {
                    share_of_slaves = floor(shares[i]) < slaves_left ?  floor(shares[i]) : slaves_left;
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
        for(int your_master = 0 ; your_master < inter_comm.size(); your_master++) {

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
        // Now that the slaves know who's their daddy, tell who are their brothers
        for(std::vector<int>::iterator slave = my_slaves.begin(); slave != my_slaves.end(); slave++) {
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
    //cout << inter_comm.rank() <<  "mumps.n " << mumps.n << endl;
    //if(inter_comm.rank() == 0){
        //cout << IRANK << " " << "mumps.n " << mumps.n << endl;
        //cout << IRANK << " " << "mumps.nz " << mumps.nz << endl;
        //exit(0);
    //}
    //if(inter_comm.rank() == 3){
        //setMumpsIcntl(4, 2);
        //setMumpsIcntl(1, 6);
        //setMumpsIcntl(2, 6);
        //setMumpsIcntl(3, 6);
    //}
    dmumps_c(&mumps);
    t = MPI_Wtime() - t;
    //cout << inter_comm.rank() << " " << partitions[0].dim(0) << " " << partitions[0].dim(1)  << " "
        //<< mumps.n << " " << mumps.nz << " " << partitions[0].NumNonzeros() << " "
        //<< mumps.write_problem << endl;
    //inter_comm.barrier();
    //exit(0);

    if(getMumpsInfo(1) != 0){
        //cout << string(32, '-') << endl
             //<< "| MUMPS ANALYSIS FAILED on MA " << setw(7) << inter_comm.rank() << " |" << endl
             //<< string(32, '-') << endl
             //<< "| info(1)       : " << setw(6) << getMumpsInfo(1) << string(4, ' ') << " |" << endl
             //<< "| info(2)       : " << setw(6) << getMumpsInfo(2) << string(4, ' ') << " |" << endl
             //<< string(32, '-') << endl;
        throw 100 * getMumpsInfo(1) - getMumpsInfo(2);
    }

    if(instance_type == 0 && verbose == true) {
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
    //if(inter_comm.rank() == 3){
        //setMumpsIcntl(4, 2);
        //setMumpsIcntl(1, 6);
        //setMumpsIcntl(2, 6);
        //setMumpsIcntl(3, 6);
    //}
    dmumps_c(&mumps);
    t = MPI_Wtime() - t;

    if(getMumpsInfo(1) != 0){
        cout << string(32, '-') << endl
             << "| MUMPS Factoriz FAILED on MA " << setw(7) << inter_comm.rank() << " |" << endl
             << string(32, '-') << endl
             << "| info(1)       : " << setw(6) << getMumpsInfo(1) << string(4, ' ') << " |" << endl
             << "| info(2)       : " << setw(6) << getMumpsInfo(2) << string(4, ' ') << " |" << endl
             << string(32, '-') << endl;
        throw 100 * getMumpsInfo(1) - getMumpsInfo(2);
    }


    if(instance_type == 0 && verbose == true) {
        double flop = getMumpsRinfoG(3);
        int prec = cout.precision();
        cout.precision(2);
        cout << string(32, '-') << endl
             << "| MUMPS FACTORIZ on MA " << setw(7) << inter_comm.rank() << " |" << endl
             << string(32, '-') << endl
             << "| N             : " << setw(12) << mumps.n << " |" << endl
             << "| NZ            : " << setw(12) << mumps.nz << " |" << endl
             << "| Flops         : " << setw(6) << scientific << flop << string(4, ' ') << " |" << endl
             << "| Time          : " << setw(6) << t << " sec |" << endl
             << "| memory        : " << setw(6) << getMumpsInfo(22) << " M|" << endl
             << string(32, '-') << endl;
        cout.precision(prec);
    }
}

MV_ColMat_double abcd::sumProject(double alpha, MV_ColMat_double &Rhs, double beta, MV_ColMat_double &X)
{
    mpi::communicator world;
    //int s = X.dim(1);
    if (alpha!=0 && beta!=0) assert(X.dim(1) == Rhs.dim(1));

    int s = alpha != 0 ? Rhs.dim(1) : X.dim(1);

    // Build the mumps rhs
    mumps.rhs = new double[mumps.n * s];
    for(int i = 0; i < mumps.n * s; i++) mumps.rhs[i] = 0;
    MV_ColMat_double mumps_rhs(mumps.rhs, mumps.n, s, MV_Matrix_::ref);

    int pos = 0;
    int b_pos = 0;
    MV_ColMat_double Delta(n, s, 0);

    if(beta != 0 || alpha != 0){
        for(int k = 0; k < partitions.size(); k++) {
            MV_ColMat_double r(partitions[k].dim(0), s, 0);
            MV_ColMat_double compressed_x(partitions[k].dim(1), s, 0);

            // avoid useless operations
            if(beta != 0){
                int x_pos = 0;
                for(int i = 0; i < local_column_index[k].size(); i++) {
                    int ci = local_column_index[k][i];
                    for(int j = 0; j < s; j++) {
                        //assert(x_pos < compressed_x.dim(0));
                        //assert(ci < X.dim(0));
                        compressed_x(x_pos, j) = X(ci, j);
                    }
                    x_pos++;
                }

                r = smv(partitions[k], compressed_x) * beta;

            }

            if(alpha != 0){
                MV_ColMat_double rr(partitions[k].dim(0), s);
                rr = Rhs(MV_VecIndex(b_pos, b_pos + partitions[k].dim(0) - 1), MV_VecIndex(0, s -1));
                r = r + rr * alpha;
            }

            b_pos += partitions[k].dim(0);

            for(int r_p = 0; r_p < s; r_p++) {
                //for(int i = pos; i < pos + partitions[k].dim(1); i++) {
                    //mumps.rhs[i + r_p * mumps.n] = 0;
                //}
                int j = 0;
                for(int i = pos + partitions[k].dim(1); i < pos + partitions[k].dim(1) + partitions[k].dim(0); i++) {
                    mumps_rhs(i, r_p) = r(j++, r_p);
                }
            }

            pos += partitions[k].dim(1) + partitions[k].dim(0);

        }

        if(infNorm(X) != 0 || infNorm(Rhs) != 0){

            int job = 1;
            mpi::broadcast(intra_comm, job, 0);

            mumps.nrhs = s;
            mumps.lrhs = mumps.n;
            mumps.job = 3;

            double t = MPI_Wtime();
            dmumps_c(&mumps);
            t = MPI_Wtime() - t;

            //cout << "[" << inter_comm.rank() << "] Time spent in direct solver : " << t << endl;

            //MV_ColMat_double Sol(mumps.rhs, mumps.n, s);

            int x_pos = 0;
            for(int k = 0; k < partitions.size(); k++) {
                for(int i = 0; i < local_column_index[k].size(); i++) {
                    int ci = local_column_index[k][i];
                    for(int j = 0; j < s; j++) {
                        Delta(ci, j) = Delta(ci, j) + mumps_rhs(x_pos, j) ;
                    }
                    x_pos++;
                }
                x_pos += partitions[k].dim(0);
            }
        }
    }

    // Where the other Deltas are going to be summed
    //MV_ColMat_double Others(n, s, 0);

    std::vector<double *>itcp(nbparts);
    std::vector<double *>otcp(nbparts);

    std::vector<mpi::status> sts;
    std::vector<mpi::request> reqs;
    int id_to = 0;

    double t = MPI_Wtime();
    double t1 = t;
    double t2 = 0;
    for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
            it != col_interconnections.end(); it++) {

        if(it->second.size() == 0) continue;
        // Prepare the data to be sent
        //
        //create it if does not exist
        itcp[it->first] = new double[it->second.size()*s];
        otcp[it->first] = new double[it->second.size()*s];

        double t3 = MPI_Wtime();
        int kp = 0;
        for(int j = 0; j < s; j++) {
            for(std::vector<int>::iterator i = it->second.begin(); i != it->second.end(); i++) {
                //itc[it->first].push_back(Delta(*i, j));
                itcp[it->first][kp] = Delta(*i,j);
                kp++;
            }
        }
        t2 += MPI_Wtime() - t3;

        reqs.push_back(inter_comm.irecv(it->first, 31, otcp[it->first], kp));
        reqs.push_back(inter_comm.isend(it->first, 31, itcp[it->first], kp));
    }

    mpi::wait_all(reqs.begin(), reqs.end());
    t1 = MPI_Wtime() - t1;

    for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
            it != col_interconnections.end(); it++) {

        if(it->second.size() == 0) continue;
        int p = 0;
        for(int j = 0; j < s; j++) {
            for(std::vector<int>::iterator i = it->second.begin(); i != it->second.end(); i++) {
                //Others(*i, j) += otcp[it->first][p];
                Delta(*i, j) += otcp[it->first][p];
                p++;
            }
        }

    }

    // Now sum the data to Delta
    //Delta += Others;

    //clog << "["<< IRANK << "] Time spent merging results : " << MPI_Wtime() -t << " ["<<t1 - t2<<", "<< t2<< "]" << endl;
    for(int i = 0; i < itcp.size(); i ++) {
        delete[] itcp[i];
        delete[] otcp[i];
    }

    delete[] mumps.rhs;

    return Delta;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::waitForSolve
 *  Description:  
 * =====================================================================================
 */
    void
abcd::waitForSolve()
{
    DMUMPS_STRUC_C mu;
    mpi::communicator world;

    int job = 0;
    do{
        mpi::broadcast(intra_comm, job, 0);
        
        if(job == -1) break;
        if(job == 1){
            mumps.job = 3;
            dmumps_c(&mumps);
        } else if (job == 2) {
            mu.sym = 2;
            mu.par = 1;
            mu.job = -1;

            mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) world);
            dmumps_c(&mu);

            mu.icntl[0] = -1;
            mu.icntl[1] = -1;
            mu.icntl[2] = -1;

            mu.job = 1;
            dmumps_c(&mu);
            mu.job = 2;
            dmumps_c(&mu);
            mu.job = 3;
            dmumps_c(&mu);
        }
    }while(true);

}		/* -----  end of function abcd::waitForSolve()  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::simpleProject
 *  Description:  Computes sumproject for only one partition without waiting
 *  for interconnections with others
 * =====================================================================================
 */
MV_ColMat_double abcd::simpleProject(MV_ColMat_double &X)
{
    int s = X.dim(1);

    // Build the mumps rhs
    mumps.rhs = new double[mumps.n * s];
    for(int i = 0; i < mumps.n * s; i++) mumps.rhs[i] = 0;
    int pos = 0;
    for(int k = 0; k < partitions.size(); k++) {
        MV_ColMat_double r(partitions[k].dim(0), s, 0);
        MV_ColMat_double compressed_x(partitions[k].dim(1), s, 0);

        int x_pos = 0;
        for(int i = 0; i < local_column_index[k].size(); i++) {
            int ci = local_column_index[k][i];
            for(int j = 0; j < s; j++) {
                //assert(x_pos < compressed_x.dim(0));
                //assert(ci < X.dim(0));
                compressed_x(x_pos, j) = X(ci, j);
            }
            x_pos++;
        }

        r = smv(partitions[k], compressed_x);
        

        for(int r_p = 0; r_p < s; r_p++) {
            for(int i = pos; i < pos + partitions[k].dim(1); i++) {
                mumps.rhs[i + r_p * mumps.n] = 0;
            }
            int j = 0;
            for(int i = pos + partitions[k].dim(1); i < pos + partitions[k].dim(1) + partitions[k].dim(0); i++) {
                mumps.rhs[i + r_p * mumps.n] = r(j++, r_p);
            }
        }

        pos += partitions[k].dim(1) + partitions[k].dim(0);

    }

    int job = 1;
    mpi::broadcast(intra_comm, job, 0);

    mumps.nrhs = s;
    mumps.lrhs = mumps.n;
    mumps.job = 3;

    dmumps_c(&mumps);

    MV_ColMat_double Sol(mumps.rhs, mumps.n, s);
    MV_ColMat_double Delta(n, s, 0);

    int x_pos = 0;
    for(int k = 0; k < partitions.size(); k++) {
        for(int i = 0; i < local_column_index[k].size(); i++) {
            int ci = local_column_index[k][i];
            for(int j = 0; j < s; j++) {
                Delta(ci, j) = Sol(x_pos, j) ;
            }
            x_pos++;
        }
        x_pos += partitions[k].dim(0);
    }
    delete[] mumps.rhs;

    return Delta;
}

MV_ColMat_double abcd::coupleSumProject(double alpha, MV_ColMat_double &Rhs, double beta, MV_ColMat_double &X, int my_bro)
{
    mpi::communicator world;
    //int s = X.dim(1);
    if (alpha!=0 && beta!=0) assert(X.dim(1) == Rhs.dim(1));

    int s = alpha != 0 ? Rhs.dim(1) : X.dim(1);

    // Build the mumps rhs
    mumps.rhs = new double[mumps.n * s];
    for(int i = 0; i < mumps.n * s; i++) mumps.rhs[i] = 0;
    int pos = 0;
    int b_pos = 0;
    MV_ColMat_double Delta(n, s, 0);

    if(beta != 0 || alpha != 0){
        for(int k = 0; k < partitions.size(); k++) {
            MV_ColMat_double r(partitions[k].dim(0), s, 0);
            MV_ColMat_double compressed_x(partitions[k].dim(1), s, 0);

            // avoid useless operations
            if(beta != 0){
                int x_pos = 0;
                for(int i = 0; i < local_column_index[k].size(); i++) {
                    int ci = local_column_index[k][i];
                    for(int j = 0; j < s; j++) {
                        //assert(x_pos < compressed_x.dim(0));
                        //assert(ci < X.dim(0));
                        compressed_x(x_pos, j) = X(ci, j);
                    }
                    x_pos++;
                }

                r = smv(partitions[k], compressed_x) * beta;
            }


            if(alpha != 0){
                MV_ColMat_double rr(partitions[k].dim(0), s);
                rr = Rhs(MV_VecIndex(b_pos, b_pos + partitions[k].dim(0) - 1), MV_VecIndex(0, s -1));
                r = r + rr * alpha;
            }
            
            b_pos += partitions[k].dim(0);

            for(int r_p = 0; r_p < s; r_p++) {
                for(int i = pos; i < pos + partitions[k].dim(1); i++) {
                    mumps.rhs[i + r_p * mumps.n] = 0;
                }
                int j = 0;
                for(int i = pos + partitions[k].dim(1); i < pos + partitions[k].dim(1) + partitions[k].dim(0); i++) {
                    mumps.rhs[i + r_p * mumps.n] = r(j++, r_p);
                }
            }

            pos += partitions[k].dim(1) + partitions[k].dim(0);

        }

        if(infNorm(X) != 0 || infNorm(Rhs) != 0){

            int job = 1;
            mpi::broadcast(intra_comm, job, 0);

            mumps.nrhs = s;
            mumps.lrhs = mumps.n;
            mumps.job = 3;

            double t = MPI_Wtime();
            dmumps_c(&mumps);
            t = MPI_Wtime() - t;

            //cout << "[" << inter_comm.rank() << "] Time spent in direct solver : "
                //<< t << endl;

            MV_ColMat_double Sol(mumps.rhs, mumps.n, s);

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
        }
    }


    // Where the other Deltas are going to be summed
    MV_ColMat_double Others(n, s, 0);


    std::map<int, std::vector<double> > itc;
    std::map<int, std::vector<double> > otc;
    std::vector<mpi::status> sts;
    mpi::request *reqs = new mpi::request[2*col_interconnections.size()];
    int id_to = 0;

    for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
            it != col_interconnections.end(); it++) {

        if(it->first != my_bro) continue;
        // Prepare the data to be sent
        //
        //create it if does not exist
        itc[it->first];
        int kp = 0;
        for(int j = 0; j < s; j++) {
            for(std::vector<int>::iterator i = it->second.begin(); i != it->second.end(); i++) {
                itc[it->first].push_back(Delta(*i, j));
                kp++;
            }
        }

        otc[it->first];
        reqs[id_to++] = inter_comm.irecv(it->first, 31, otc[it->first]);
        reqs[id_to++] = inter_comm.isend(it->first, 31, itc[it->first]);
    }

    mpi::wait_all(reqs, reqs + id_to);
    delete[] reqs;
    for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
            it != col_interconnections.end(); it++) {

        if(it->first != my_bro) continue;
        int p = 0;
        for(int j = 0; j < s; j++) {
            for(std::vector<int>::iterator i = it->second.begin(); i != it->second.end(); i++) {
                Others(*i, j) += otc[it->first][p];
                p++;
            }
        }

    }
    // Now sum the data to Delta
    Delta += Others;
    delete[] mumps.rhs;

    return Delta;
}

MV_ColMat_double abcd::spSimpleProject(std::vector<int> mycols)
{
    bool dense_rhs = (icntl[13] == 1);
    //dense_rhs = true;

    int s = mycols.size();
    // Build the mumps rhs

    mumps.rhs = new double[mumps.n * s];
    for(int i = 0; i < mumps.n * s; i++) mumps.rhs[i] = 0;

    MV_ColMat_double mumps_rhs(mumps.rhs, mumps.n, s, MV_Matrix_::ref);

    CompCol_Mat_double mumps_comp_rhs;

    std::vector<int> rr;
    std::vector<double> rv;

    std::vector<CompRow_Mat_double> r;
    std::vector<std::map<int,int> > loc_cols(partitions.size());

    int nzr_estim = 0;

    for(int k = 0; k < partitions.size(); k++) {

        CompRow_Mat_double Y;

        VECTOR_int yr(mycols.size(), 0);
        VECTOR_int yc(mycols.size(), 0);
        VECTOR_double yv(mycols.size(), 0);
        int c;

        int ct = 0;
        for(int i = 0; i < mycols.size(); i++){
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

        nzr_estim += Y.NumNonzeros();

        r.push_back( spmm(partitions[k], Y) );
    }


    if(dense_rhs){
        // dense mumps rhs
        mumps.icntl[20 - 1] = 0;

        int pos = 0;
        for(int k = 0; k < partitions.size(); k++) {
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
        mumps.icntl[20 - 1] = 1;

        {
            mumps.irhs_ptr      = new int[s + 1];
            rr.reserve(nzr_estim);
            rv.reserve(nzr_estim);

            int cnz = 1;
            for(int r_p = 0; r_p < s; r_p++) {
                mumps.irhs_ptr[r_p] = cnz;
                int pos = 0;
                for(int k = 0; k < partitions.size(); k++) {
                    // no need to put zeros before!
                    //
                    int j = 0;
                    for(int i = pos + partitions[k].dim(1); i < pos + partitions[k].dim(1) + partitions[k].dim(0); i++) {
                        if(r[k](j, r_p) != 0){
                            rr.push_back(i+1);
                            rv.push_back(r[k](j, r_p));
                            cnz++;
                        }
                        j++;
                    }
                    pos += partitions[k].dim(1) + partitions[k].dim(0);
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

    MV_ColMat_double Delta(size_c, s, 0);

    int x_pos = 0;

    for(int k = 0; k < partitions.size(); k++) {
        int start_c = glob_to_part[k][stC[k]];
        x_pos += glob_to_part[k][stC[k]];

        for(int i = start_c; i < column_index[k].size(); i++){

            int ci = column_index[k][i] - n_o;

            for(int j = 0; j < s; j++) {
                Delta(ci, j) = Delta(ci, j) - mumps_rhs(x_pos, j) ; // Delta = - \sum (sol)
            }

            x_pos++;
        }

        for(int j = 0; j < s; j++) {
            if(loc_cols[k][j]){
                int c = mycols[j];

                Delta(c, j) = 0.5 + Delta(c,j);
            }
        }

        x_pos += partitions[k].dim(0);
    }

    // disable sparse mumps rhs
    mumps.icntl[20 - 1] = 0;

    //delete mumps.rhs;
    delete[] mumps.rhs;

    return Delta;
}
