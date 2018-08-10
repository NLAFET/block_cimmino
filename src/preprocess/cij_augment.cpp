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
 * \file cij_augment.cpp
 * \brief Implementation of the augmentation with blocks Cij/-I
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include "abcd.h"

/*!
 *  \brief Augments the matrix using Cij/-I blocks
 *
 *  Augment the matrix M by building the augmentation part C with Cij/-I blocks
 *  on the interconnected columns. If filtering of the augmentation is triggered
 *  this function also filters columns.
 *
 *  \param M: A vector of the compressed partitions (splib format) to augment
 *
 */
void abcd::cijAugmentMatrix(std::vector<CompCol_Mat_double> &M)
{
    int nbcols = A.dim(1); // size of Abar
    int nz_c = 0;
    std::map<int,std::vector<CompCol_Mat_double> > C; // Augmentation blocks
    /* stCols(i) contains the number of augmented columns per other interconnected partition */
    std::map<int,std::vector<int> > stCols;
#ifdef WIP
    double filter_c = dcntl[Controls::aug_filter];
#else
    double filter_c = 0;
#endif //WIP
    stC = std::vector<int>(M.size(), -1);

    if (icntl[Controls::scaling] == 0)
        LWARNING << "Using C_ij based augmentation without scaling could induce numerical instability";

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

            // Build Cij/-I submatrices from interconnections
            CompCol_Mat_double C_ij;
            {
                CompRow_Mat_double A_ij(sub_matrix(M[i], intersect));
                CompRow_Mat_double A_ji(sub_matrix(M[j], intersect));
                CompRow_Mat_double A_jiT = csr_transpose(A_ji);
                C_ij = spmm(A_ij, A_jiT);
            }

            if(C_ij.NumNonzeros() == 0) continue;

            /* Compress C_ij */
            std::vector<int> cic, cir, ci;
            bool reversed = false;
            std::vector<int> selected_cols;

            // Compute the Column compression
            int l = 0;
            for(int k = 1; k <= C_ij.dim(1); k++) {
                if(C_ij.col_ptr(k) != C_ij.col_ptr(k - 1)){
                    bool valid = false;
                    int coli = C_ij.col_ptr(k-1);
                    while(coli < C_ij.col_ptr(k)){
#ifdef WIP
                        // Select columns of S: keep if value exists superior to threshold ????
                        if(icntl[Controls::aug_iterative] != 0){ 
                            if(abs(C_ij.val(coli)) >= dcntl[Controls::aug_precond]){
                                selected_S_columns.push_back( nbcols + k - n_o);
                            } else {
                                skipped_S_columns.push_back( nbcols + k - n_o);
                            }
                        }
#endif //WIP

                        // If value exists superior to threshold, keep the column ???
                        if(abs(C_ij.val(coli)) >= filter_c){
                            valid = true;
                            break;
                        }

#ifdef WIP
                        // If we solve S iterative, then no filtering needed, only selected columns
                        if( icntl[Controls::aug_iterative] != 2 ){ // don't reduce, we just need the selected columns!
                            valid = true; // let the force be with you, always!
                            break;
                        }
#endif //WIP
                        coli++;
                    }
                    if(valid) cic.push_back(l);
                }
                l++;
            }

            // Compute the Row compression
            CompCol_Mat_double CT_ij = csc_transpose(C_ij);
            l = 0;
            for(int k = 1; k <= CT_ij.dim(1); k++) {
                if(CT_ij.col_ptr(k) != CT_ij.col_ptr(k - 1)){
                    bool valid = false;
                    int coli = CT_ij.col_ptr(k-1);
                    while(coli < CT_ij.col_ptr(k)){
#ifdef WIP
                        // Select rows of S: keep if value exists superior to threshold ????
                        if(icntl[Controls::aug_iterative] != 0){
                            if(abs(CT_ij.val(coli)) >= dcntl[Controls::aug_precond]){
                                selected_S_columns.push_back( nbcols + k - n_o);
                            } else {
                                skipped_S_columns.push_back( nbcols + k - n_o);
                            }
                        }

#endif //WIP
                        // If value exists superior to threshold, keep the row ???
                        if(abs(CT_ij.val(coli)) >= filter_c){
                            valid = true;
                            break;
                        }

#ifdef WIP
                        // If we solve S iterative, then no filtering needed, only selected columns
                        if( icntl[Controls::aug_iterative] != 2 ){ // don't reduce, we just need the selected columns!
                            valid = true; // let the force be with you, always!
                            break;
                        }
#endif //WIP
                        coli++;
                    }
                    if(valid) cir.push_back(l);
                }
                l++;
            }


            // If we have less rows than columns, then transpose C_ij
            if(cic.size() <= cir.size()) {
                ci = cic;
            } else {
                ci = cir;
                C_ij = CT_ij;
                reversed = true;
            }

            if(ci.empty()) continue;

            int n_cij_before = C_ij.dim(1);

            // Finally get the filtered Augmentation blocks
            C_ij = sub_matrix(C_ij, ci);

            //// Build compressed I
            CompCol_Mat_double I;
            {
                VECTOR_int ir(C_ij.dim(1));
                VECTOR_int ic(C_ij.dim(1)+1);
                VECTOR_double iv(C_ij.dim(1));
                for(int k = 0; k < C_ij.dim(1); k++){
                    ir(k) = ci[k];
                    ic(k) = k;
                    iv(k) = -1.0;
                }
                ic(C_ij.dim(1)) = C_ij.dim(1) + 1;

                I = CompCol_Mat_double(n_cij_before, C_ij.dim(1), C_ij.dim(1),
                        iv, ir, ic);
            }

            stCols[i].push_back(nbcols);
            stCols[j].push_back(nbcols);

            if(!reversed){
                C[i].push_back(C_ij);
                C[j].push_back(I);
            } else {
                C[i].push_back(I);
                C[j].push_back(C_ij);
            }

            nz_c += C_ij.NumNonzeros() + I.NumNonzeros();
            nbcols += C_ij.dim(1);
        }

    }
    // Size of C
    size_c = nbcols - A.dim(1);
    n = nbcols;
    LINFO << "Size of C : " << size_c;

#ifdef WIP
    if(icntl[Controls::aug_analysis] != 0) return;
#endif //WIP

    // Augment the matrices
    for(size_t k = 0; k < M.size(); k++){
        if(stCols[k].size() == 0) continue;
        // now augment each partition!
        stC[k] = stCols[k][0];
        M[k] = concat_columns(M[k], C[k], stCols[k]);
        M[k] = resize_columns(M[k], nbcols);
    }
}               /* -----  end of function abcd::cijAugmentMatrix  ----- */
