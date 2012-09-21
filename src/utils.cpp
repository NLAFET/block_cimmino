#include "abcd.h"

/// Partition weigts
void abcd::partitionWeights(std::vector<int> &parts, std::vector<int> weights, int nb_parts)
{
    int mean = std::accumulate(weights.begin(), weights.end(), 0) / nb_parts;
    int cum = 0;
    int precum = 0;

    for(int c = 0; c < weights.size(); c++) {
        precum = cum;
        cum += weights[c];

        if(cum >= mean) {
            if((mean - precum) > (cum - mean)) {
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
Eigen::MatrixXd abcd::ddot(Eigen::MatrixXd p, Eigen::MatrixXd ap)
{
    int ln = p.cols();
    int lm = p.rows();
    int rn = ap.cols();
    int rm = ap.rows();
    if(lm != rm) throw - 800;

    Eigen::MatrixXd loc_p(lm, ln);
    Eigen::MatrixXd loc_ap(rm, rn);
    Eigen::MatrixXd loc_r(ln, rn);
    Eigen::MatrixXd r(ln, rn);
    loc_ap.setZero();
    loc_p.setZero();
    loc_r.setZero();

    int pos = 0;
    for(int i = 0; i < lm; i++) {
        if(comm_map(i) == 1) {
            for(int j = 0; j < ln; j++) {
                loc_p(pos, j) = p(i, j);
            }
            for(int j = 0; j < rn; j++) {
                loc_ap(pos, j) = ap(i, j);
            }
            pos++;
        }
    }

    // R = P'AP
    loc_r = loc_p.transpose() * loc_ap;

    const double *l_r_ptr = loc_r.data();
    double *r_ptr = r.data();
    mpi::all_reduce(inter_comm, l_r_ptr, ln * rn,  r_ptr, std::plus<double>());

    return r;
}

void abcd::get_nrmres(Eigen::MatrixXd v, double &nrmR, double &nrmX)
{
    int rn = v.cols();
    int rm = v.rows();
    nrmX = 0;
    nrmR = 0;

    Eigen::MatrixXd loc_v(rm, rn);
    Eigen::MatrixXd loc_r(m, rn);
    loc_v.setZero();
    loc_r.setZero();

    int pos = 0;
    for(int i = 0; i < rm; i++) {
        if(comm_map(i) == 1) {
            for(int j = 0; j < rn; j++) {
                //@todo : this stands only when there is one rhs
                nrmX += pow(v(i, j), 2);
            }
            pos++;
        }
    }

    pos = 0;
    for(int p = 0; p < parts.size(); p++) {
        Eigen::VectorXd compressed_x(parts[p].cols());
        for(int j = 0; j < v.cols(); j++) {
            compressed_x.setZero();

            int x_pos = 0;
            for(int i = 0; i < local_column_index[p].size(); i++) {
                int ci = local_column_index[p][i];
                compressed_x(x_pos) = v(ci, j);
                x_pos++;
            }
            loc_r.col(j).segment(pos, parts[p].rows()) = parts[p] * compressed_x;
        }
        //nrmX += compressed_x.squaredNorm();

        pos += parts[p].rows();
    }

    loc_r  = b - loc_r;
    //nrmX = v.squaredNorm();

    double loc_nrm = loc_r.squaredNorm();
    //double loc_nrm = loc_r.lpNorm<Infinity>();
    double nrm;
    mpi::all_reduce(inter_comm, &loc_nrm, 1,  &nrm, std::plus<double>());
    nrmR = sqrt(nrm);
    mpi::all_reduce(inter_comm, &nrmX, 1,  &nrm, std::plus<double>());
    nrmX = sqrt(nrm);
}

/// Compair pairs
bool ip_comp(const dipair &l, const dipair &r)
{
    return l.first > r.first;
}

/// Sum nonzeros from parts
int sum_nnz(int res, Eigen::SparseMatrix<double, RowMajor> M)
{
    return res += M.nonZeros();
}
int sum_rows(int res, Eigen::SparseMatrix<double, RowMajor> M)
{
    return res += M.rows();
}
int sum_cols(int res, Eigen::SparseMatrix<double, RowMajor> M)
{
    return res += M.cols();
}
bool comp_cols(Eigen::SparseMatrix<double, RowMajor> L, Eigen::SparseMatrix<double, RowMajor> R)
{
    return (L.cols() < R.cols());
}
