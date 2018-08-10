// Copyright Institut National Polytechnique de Toulouse (2014) 
// Contributor(s) :
// M. Zenadi <mzenadi@enseeiht.fr>
// D. Ruiz <ruiz@enseeiht.fr>
// R. Guivarch <guivarch@enseeiht.fr>

// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html"

// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 

// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 

// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.

/*!
 * \file aij_augment.cpp
 * \brief Implementation of the augmentation with blocks Aij/-Aij
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include "abcd.h"

/*!
 *  \brief Augments the matrix using Aij/-Aji blocks
 *
 *  Augment the matrix M by building the augmentation part C with Aij/-Aji blocks
 *  on the interconnected columns. If filtering of the augmentation is triggered
 *  this function also filters columns.
 *
 *  \param M: A vector of the compressed partitions (splib format) to augment
 *
 */
void abcd::aijAugmentMatrix(std::vector<CompCol_Mat_double> &M)
{

    int nbcols = A.dim(1); // size of Abar
    std::map<int,std::vector<CompCol_Mat_double> > C; // Augmentation blocks
    /* stCols(i) contains the number of augmented columns per other interconnected partition */
    std::map<int,std::vector<int> > stCols;
    stC.assign(M.size(), -1);

    // For each couple of partitions build the augmentation
    for( size_t i = 0; i < M.size() - 1; i++ ){
        for ( size_t j = i+1; j < M.size(); j++ ) {
            // find interconnected columns between partitions
            std::vector<int> intersect;
            std::set_intersection(column_index[i].begin(),
                                  column_index[i].end(),
                                  column_index[j].begin(),
                                  column_index[j].end(),
                                  std::back_inserter(intersect));

            if (intersect.empty()) continue;

            // Build Aij/-Aji submatrices from interconnections
            CompCol_Mat_double A_ij(sub_matrix(M[i], intersect));
            CompCol_Mat_double A_ji(sub_matrix(M[j], intersect));
            double *jv = A_ji.val_ptr();
            for (int k = 0; k < A_ji.NumNonzeros(); k++) {
                jv[k] *= -1.0;
            }

#ifdef WIP
            if(dcntl[Controls::aug_filter] != 0 || icntl[Controls::aug_iterative] != 0) {
                std::vector<int> selected_cols;
                std::vector<double> frob_ij, mu;

                frob_ij.reserve(A_ij.dim(1));
                mu.reserve(A_ij.dim(1));
                double card_max = 0;
                double frob_sum = 0;
                double nu;

                for (int k = 0; k < A_ij.dim(1); ++k){
                    // Column k of the augmentation
                    VECTOR_int A_ij_k_ind, A_ji_k_ind;
                    VECTOR_double A_ij_k = middleCol(A_ij, k, A_ij_k_ind);
                    VECTOR_double A_ji_k = middleCol(A_ji, k, A_ji_k_ind);

                    double card_current = A_ij_k_ind.size() * A_ji_k_ind.size();

                    // Compute the Froebenius norm of each rank one update corresponding to column k
                    // exploit the sparcity of the vectors!
                    frob_ij.push_back(sqrt( squaredNorm(A_ij_k, A_ij_k_ind) * squaredNorm(A_ji_k, A_ji_k_ind)));
                    frob_sum += frob_ij[k];

                    // Underestimate of the number of entries in AijAji^T
                    card_max = card_max > card_current ? card_max : card_current;
                }

                // Mean distribution of the rank one contributions
                nu = (frob_sum / frob_ij.size()) / sqrt(card_max);

                for (int k = 0; k < A_ij.dim(1); ++k){
                    // Column k of the augmentation
                    VECTOR_int A_ij_k_ind, A_ji_k_ind;
                    VECTOR_double A_ij_k = middleCol(A_ij, k, A_ij_k_ind);
                    VECTOR_double A_ji_k = middleCol(A_ji, k, A_ji_k_ind);

                    // Count the number of influential entries
                    double inf_ij = infNorm(A_ij_k);
                    double inf_ji = infNorm(A_ji_k);

                    double p = 0, q = 0; // card_ij, card_ji
                    for (int l = 0; l < A_ij_k_ind.size(); ++l){
                        if (abs(A_ij_k(A_ij_k_ind(l))) >= nu/inf_ji) p++;
                    }
                    for (int l = 0; l < A_ji_k_ind.size(); ++l){
                        if (abs(A_ji_k(A_ji_k_ind(l))) >= nu/inf_ij) q++;
                    }

                    // In case the cardinality is zero, we just set to #rows in Aij/Aji
                    p = ( p==0 ? A_ij_k_ind.size() : p );
                    q = ( q==0 ? A_ji_k_ind.size() : q );

                    // Define a scaled norm
                    double mu_ij_k = frob_ij[k] / sqrt(p*q);
                    mu.push_back(mu_ij_k);

                    // Select columns depending on threshold
                    if(mu_ij_k >= dcntl[Controls::aug_filter]){
                        selected_cols.push_back(k);
                    }

                    // Select columns of S: keep if value exists superior to threshold ????
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

                // Finally get the filtered Augmentation blocks
                if( icntl[Controls::aug_iterative] != 2 ) { // don't reduce the A_ij/A_ji, we just need the selected columns!
                    A_ij = sub_matrix(A_ij, selected_cols);
                    A_ji = sub_matrix(A_ji, selected_cols);
                }
            }
#endif //WIP

            stCols[i].push_back(nbcols);
            stCols[j].push_back(nbcols);
            C[i].push_back(A_ij);
            C[j].push_back(A_ji);

            nbcols += A_ij.dim(1);
        }
    }

    // Size of C
    size_c = nbcols - A.dim(1);
    n = nbcols;
    LINFO << "Size of C : " << size_c;

#ifdef WIP
    if(icntl[Controls::aug_analysis] != 0) return;
#endif // WIP

    // Augment the matrices
    for(size_t k = 0; k < M.size(); k++){
        if(stCols[k].size() == 0) continue;
        // now augment each partition!
        stC[k] = stCols[k][0];
        M[k] = concat_columns(M[k], C[k], stCols[k]);
        M[k] = resize_columns(M[k], nbcols);
    }
}               /* -----  end of function abcd::aijAugmentMatrix  ----- */
