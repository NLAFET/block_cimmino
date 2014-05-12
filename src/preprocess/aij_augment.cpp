#include "abcd.h"

void abcd::aijAugmentMatrix(std::vector<CompCol_Mat_double> &M)
{

    int nbcols = A.dim(1);
    std::map<int,std::vector<CompCol_Mat_double> > C;
    std::map<int,std::vector<int> > stCols;
    double filter_c = dcntl[Controls::aug_filter];
    stC = std::vector<int>(M.size(), -1);

    for( size_t i = 0; i < M.size() - 1; i++ ){
        for ( size_t j = i+1; j < M.size(); j++ ) {
            std::vector<int> intersect;
            std::set_intersection(column_index[i].begin(), column_index[i].end(),
                                    column_index[j].begin(), column_index[j].end(),
                                    std::back_inserter(intersect));
            if (intersect.empty()) continue;

            CompCol_Mat_double A_ij = sub_matrix(M[i], intersect);
            CompCol_Mat_double A_ji = sub_matrix(M[j], intersect);

            for(int k = 0; k < A_ji.NumNonzeros(); k++)
                A_ji.val(k) *= double(-1);

            if(filter_c != 0 || icntl[Controls::aug_iterative] != 0) {
                std::vector<int> selected_cols;
                std::vector<double> frob_ij, mu;

                frob_ij.reserve(A_ij.dim(1));
                mu.reserve(A_ij.dim(1));
                double card_max = 0;
                double frob_sum = 0;
                double nu;

                for (int k = 0; k < A_ij.dim(1); ++k){
                    VECTOR_int A_ij_k_ind, A_ji_k_ind;
                    VECTOR_double A_ij_k = middleCol(A_ij, k, A_ij_k_ind);
                    VECTOR_double A_ji_k = middleCol(A_ji, k, A_ji_k_ind);

                    double card_current = A_ij_k_ind.size() * A_ji_k_ind.size();

                    // exploit the sparcity of the vectors!
                    frob_ij.push_back(sqrt( squaredNorm(A_ij_k, A_ij_k_ind) * squaredNorm(A_ji_k, A_ji_k_ind)));
                    frob_sum += frob_ij[k];

                    card_max = card_max > card_current ? card_max : card_current;
                }

                nu = (frob_sum / frob_ij.size()) / sqrt(card_max);

                for (int k = 0; k < A_ij.dim(1); ++k){
                    VECTOR_int A_ij_k_ind, A_ji_k_ind;
                    VECTOR_double A_ij_k = middleCol(A_ij, k, A_ij_k_ind);
                    VECTOR_double A_ji_k = middleCol(A_ji, k, A_ji_k_ind);

                    double inf_ij = infNorm(A_ij_k);
                    double inf_ji = infNorm(A_ji_k);

                    double p = 0, q = 0;
                    for (int l = 0; l < A_ij_k_ind.size(); ++l){
                        if (abs(A_ij_k(A_ij_k_ind(l))) >= nu/inf_ji) p++;
                    }
                    for (int l = 0; l < A_ji_k_ind.size(); ++l){
                        if (abs(A_ji_k(A_ji_k_ind(l))) >= nu/inf_ij) q++;
                    }

                    p = ( p==0 ? A_ij_k_ind.size() : p );
                    q = ( q==0 ? A_ji_k_ind.size() : q );

                    double mu_ij_k = frob_ij[k] / sqrt(p*q);
                    mu.push_back(mu_ij_k);

                    if(mu_ij_k >= filter_c){
                        selected_cols.push_back(k);
                    }

                    if(icntl[Controls::aug_iterative] != 0){ 
                        if (dcntl[Controls::aug_precond] < 0) {
                            if ((nbcols + k - n_o) % abs((int)dcntl[Controls::aug_precond]) == 0)
                                selected_S_columns.push_back( nbcols + k - n_o);
                            else
                                skipped_S_columns.push_back( nbcols + k - n_o);
                        } else {
                            if (mu_ij_k >= dcntl[Controls::aug_precond])
                                selected_S_columns.push_back( nbcols + k - n_o);
                            else
                                skipped_S_columns.push_back( nbcols + k - n_o);
                        }
                    }

                }

                if (selected_cols.empty()) continue;

                if( icntl[Controls::aug_iterative] != 2 ) { // don't reduce the A_ij/A_ji, we just need the selected columns!
                    A_ij = sub_matrix(A_ij, selected_cols);
                    A_ji = sub_matrix(A_ji, selected_cols);
                }
            }

            stCols[i].push_back(nbcols);
            stCols[j].push_back(nbcols);
            C[i].push_back(A_ij);
            C[j].push_back(A_ji);

            nbcols += A_ij.dim(1);
        }
    }

    size_c = nbcols - A.dim(1);
    n = nbcols;
    LINFO << "Size of C : " << size_c;

    if(icntl[Controls::aug_analysis] != 0) return;

    // Augment the matrices
    for(size_t k = 0; k < M.size(); k++){
        if(stCols[k].size() == 0) continue;
        // now augment each partition!
        stC[k] = stCols[k][0];
        M[k] = concat_columns(M[k], C[k], stCols[k]);
        M[k] = resize_columns(M[k], nbcols);
    }
}
