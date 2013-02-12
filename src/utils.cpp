#include <abcd.h>
#include <Eigen/src/misc/blas.h>

/// Partition weigts
void abcd::partitionWeights(std::vector<int> &parts, std::vector<int> weights, int nb_parts)
{
    int total_size = std::accumulate(weights.begin(), weights.end(), 0);
    int mean = total_size / nb_parts;
    int cum = 0;
    int precum = 0;

    for(int c = 0; c < weights.size(); c++) {
        precum = cum;
        cum += weights[c];

        if(cum > mean) {
            if((mean - precum) > 1.5*(cum - mean)) {
                parts.push_back(c);
                cum = 0;
            } else {
                parts.push_back(c - 1);
                cum = weights[c];
            }
        }
    }
    // the last partitions were not taken as they are less than the mean
    if(parts.size() == nb_parts) {
        if(parts[nb_parts - 1 ] < weights.size())
            parts[nb_parts - 1 ] = weights.size() - 1;
    } else {
        parts.push_back(weights.size() - 1);
    }
}
///DDOT
double abcd::ddot(VECTOR_double &p, VECTOR_double &ap)
{
    int lm = p.size();
    int rm = ap.size();
    if(lm != rm) throw - 800;

    VECTOR_double loc_p(lm, 0);
    VECTOR_double loc_ap(rm, 0);
    double loc_r, r;

    int pos = 0;
    for(int i = 0; i < lm; i++) {
        if(comm_map[i] == 1) {
            loc_p(pos) = p(i);
            loc_ap(pos) = ap(i);
            pos++;
        }
    }

    // R = P'AP
    for(int i = 0; i < lm; i++){
        loc_r += loc_p(i) * loc_ap(i);
    }

    mpi::all_reduce(inter_comm, loc_r, r, std::plus<double>());

    return r;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::get_nrmres
 *  Description:  Computes ||X_k|| and ||X_f - X_k||/||X_f||
 * =====================================================================================
 */
void abcd::get_nrmres(MV_ColMat_double &x, double &nrmR, double &nrmX, double &nrmXfmX)
{
    mpi::communicator world;
    int rn = x.dim(1);
    int rm = x.dim(0);
    nrmX = 0;
    nrmR = 0;

    MV_ColMat_double loc_x(rm, rn, 0);
    MV_ColMat_double loc_r(m, rn, 0);
    MV_ColMat_double loc_xfmx(rm, rn, 0);

    int pos = 0;
    for(int i = 0; i < rm; i++) {
        if(comm_map[i] == 1) {
            for(int j = 0; j < rn; j++) {
                //@todo : this stands only when there is one rhs
                nrmX += pow(x(i, j), 2);
            }
            pos++;
        }
    }

    pos = 0;
    for(int p = 0; p < partitions.size(); p++) {
        for(int j = 0; j < x.dim(1); j++) {
            VECTOR_double compressed_x = VECTOR_double((partitions[p].dim(1)), 0);

            int x_pos = 0;
            for(int i = 0; i < local_column_index[p].size(); i++) {
                int ci = local_column_index[p][i];
                compressed_x(x_pos) = x(ci, j);
                x_pos++;
            }
            VECTOR_double vj = loc_r(j);
            vj(MV_VecIndex(pos, pos+partitions[p].dim(0) - 1)) = partitions[p] * compressed_x;
            loc_r.setCol(vj, j);
        }

        pos += partitions[p].dim(0);
    }

    double loc_nrmxfmx;
    if(use_xf){
        if(inter_comm.size()>1) {
            if(inter_comm.rank()==0)
                cout << "NOT IMPLEMENTED YET : Parallel use of xf" << endl;
        } else {
            loc_xfmx = Xf - x;
            loc_nrmxfmx = infNorm(loc_xfmx);
        }
    }

    loc_r  = B - loc_r;

    double loc_nrm = loc_r.squaredSum();

    double nrms[2] = {loc_nrm, nrmX};
    double nrms_out[2];
    mpi::all_reduce(inter_comm, nrms, 2, nrms_out, std::plus<double>());

    nrmR = sqrt(nrms_out[0]);
    nrmX = sqrt(nrms_out[1]);

    if(use_xf){
        if(inter_comm.size()==1)
            mpi::all_reduce(inter_comm, &loc_nrmxfmx, 1,  &nrmXfmX, mpi::maximum<double>());
        //nrmXfmX = loc_xfmx.norm();
    }
}

/// Compair pairs
bool ip_comp(const dipair &l, const dipair &r)
{
    return l.first > r.first;
}


