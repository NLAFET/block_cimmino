#include<abcd.h>
#include<mumps.h>

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
    double ti = 0, to;
    double *xpt = X.ptr();
    int xlda = X.lda();
    int dlda = Delta.lda();

    if(beta != 0 || alpha != 0){
        for(int k = 0; k < nb_local_parts; k++) {

            CompRow_Mat_double *part = &partitions[k];
            std::vector<int> *loc_col_idx  = &local_column_index[k];
            int *fast_col_idx = fast_local_column_index[k];

            MV_ColMat_double r(part->dim(0), s, 0);
            double *rpt = r.ptr();
            int rlda = r.lda();
            MV_ColMat_double compressed_x(part->dim(1), s, 0);
            double *cxpt = compressed_x.ptr();
            int cxlda = compressed_x.lda();

            // avoid useless operations
            if(beta != 0){
                int x_pos = 0;
                for(int i = 0; i < loc_col_idx->size(); i++) {
                    int ci = fast_col_idx[i];
                    for(int j = 0; j < s; j++) {
                        cxpt[x_pos + j * cxlda] = xpt[ci + j * xlda];
                    }
                    x_pos++;
                }

                r = smv(*part, compressed_x) * beta;
            }

            if(alpha != 0 && beta !=0){
                MV_ColMat_double rr(part->dim(0), s);
                rr = Rhs(MV_VecIndex(b_pos, b_pos + part->dim(0) - 1), MV_VecIndex(0, s -1));
                r = r + rr * alpha;
            }

            if(alpha != 0 && beta ==0){
                r = Rhs(MV_VecIndex(b_pos, b_pos + part->dim(0) - 1), MV_VecIndex(0, s -1)) * alpha;
            }

            b_pos += part->dim(0);

            //to = MPI_Wtime();

            int j = 0;
            for(int i = pos + part->dim(1); i < pos + part->dim(1) + part->dim(0); i++) {
                //for(int r_p = 0; r_p < s; r_p++) mumps_rhs(i, r_p) = r(j, r_p);
                for(int r_p = 0; r_p < s; r_p++) mumps.rhs[i + r_p * mumps.n] = rpt[j + r_p * rlda];
                j++;
            }
            //ti += MPI_Wtime() - to;

            pos += part->dim(1) + part->dim(0);

        }

        int job = 1;
        mpi::broadcast(intra_comm, job, 0);

        mumps.nrhs = s;
        mumps.lrhs = mumps.n;
        mumps.job = 3;

        double t = MPI_Wtime();
        dmumps_c(&mumps);
        t = MPI_Wtime() - t;

        int x_pos = 0;
        double *dpt = Delta.ptr();
        if(nb_local_parts > 1)
        {
            for(int k = 0; k < nb_local_parts; k++) {
                for(int i = 0; i < local_column_index[k].size(); i++) {
                    //int ci = local_column_index[k][i];
                    int ci = fast_local_column_index[k][i];
                    for(int j = 0; j < s; j++) {
                        //Delta(ci, j) = Delta(ci, j) + mumps_rhs(x_pos, j) ;
                        dpt[ci + j * dlda] += mumps.rhs[x_pos + j * mumps.n];
                    }
                    x_pos++;
                }
                x_pos += partitions[k].dim(0);
            }

        } else
        {
            for(int i = 0; i < local_column_index[0].size(); i++) {
                int ci = fast_local_column_index[0][i];
                for(int j = 0; j < s; j++) {
                    dpt[ci + j * dlda] = mumps.rhs[x_pos + j * mumps.n];
                }
                x_pos++;
            }
        }
    }

    if(inter_comm.size() == 1) {
        delete[] mumps.rhs;
        return Delta;
    }

    // Where the other Deltas are going to be summed

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

    for(int i = 0; i < itcp.size(); i ++) {
        delete[] itcp[i];
        delete[] otcp[i];
    }

    delete[] mumps.rhs;

    return Delta;
}
