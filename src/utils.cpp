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

/// Compair pairs
bool ip_comp ( const dipair& l, const dipair& r)
   { return l.first > r.first; }
   
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
bool comp_cols(Eigen::SparseMatrix<double, RowMajor> L, Eigen::SparseMatrix<double, RowMajor> R){ return (L.cols() < R.cols()); }