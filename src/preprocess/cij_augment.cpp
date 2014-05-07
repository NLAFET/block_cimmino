#include "abcd.h"

void abcd::cijAugmentMatrix(std::vector<CompCol_Mat_double> &M)
{
    int nbcols = A.dim(1);
    int nz_c = 0;
    std::map<int,std::vector<CompCol_Mat_double> > C;
    std::map<int,std::vector<int> > stCols;
    double filter_c = dcntl[Controls::aug_filter];
    stC = std::vector<int>(M.size(), -1);

    if (icntl[Controls::scaling] == 0)
      clog << "Using C_ij based augmentation gives better results with scaling" << endl;


    for( size_t i = 0; i < M.size() - 1; i++ ){
        for ( size_t j = i+1; j < M.size(); j++ ) {
            std::vector<int> intersect;
            std::set_intersection(column_index[i].begin(), column_index[i].end(),
                                    column_index[j].begin(), column_index[j].end(),
                                    std::back_inserter(intersect));
            if (intersect.empty()) continue;

            CompCol_Mat_double C_ij;
            {
                CompRow_Mat_double A_ij = CompRow_Mat_double(sub_matrix(M[i], intersect));
                CompRow_Mat_double A_ji = CompRow_Mat_double(sub_matrix(M[j], intersect));
                CompRow_Mat_double A_jiT = csr_transpose(A_ji);
                C_ij = spmm(A_ij, A_jiT);
            }

            if(C_ij.NumNonzeros() == 0) continue;

            /*-----------------------------------------------------------------------------
             *  Compress C_ij
             *-----------------------------------------------------------------------------*/
            std::vector<int> cic, cir, ci;
            bool reversed = false;
            std::vector<int> selected_cols;

            //[> Compute the Column compression <]
            int l = 0;
            for(int k = 1; k <= C_ij.dim(1); k++) {
                if(C_ij.col_ptr(k) != C_ij.col_ptr(k - 1)){
                    bool valid = false;
                    int coli = C_ij.col_ptr(k-1);
                    while(coli < C_ij.col_ptr(k)){
                        if(icntl[Controls::aug_iterative] != 0){ 
                            if(abs(C_ij.val(coli)) >= dcntl[Controls::aug_precond]){
                                selected_S_columns.push_back( nbcols + k - n_o);
                            } else {
                                skipped_S_columns.push_back( nbcols + k - n_o);
                            }
                        }

                        if(abs(C_ij.val(coli)) >= filter_c){
                            valid = true;
                            break;
                        }

                        if( icntl[Controls::aug_iterative] != 2 ){ // don't reduce, we just need the selected columns!
                            valid = true; // let the force be with you, always!
                            break;
                        }

                        coli++;
                    }
                    if(valid) cic.push_back(l);
                }
                l++;
            }

            //[> Compute the Row compression <]
            CompCol_Mat_double CT_ij = csc_transpose(C_ij);
            l = 0;
            for(int k = 1; k <= CT_ij.dim(1); k++) {
                if(CT_ij.col_ptr(k) != CT_ij.col_ptr(k - 1)){
                    bool valid = false;
                    int coli = CT_ij.col_ptr(k-1);
                    while(coli < CT_ij.col_ptr(k)){
                        if(icntl[Controls::aug_iterative] != 0){ 
                            if(abs(CT_ij.val(coli)) >= dcntl[Controls::aug_precond]){
                                selected_S_columns.push_back( nbcols + k - n_o);
                            } else {
                                skipped_S_columns.push_back( nbcols + k - n_o);
                            }
                        }

                        if(abs(CT_ij.val(coli)) >= filter_c){
                            valid = true;
                            break;
                        }

                        if( icntl[Controls::aug_iterative] != 2 ){ // don't reduce, we just need the selected columns!
                            valid = true; // let the force be with you, always!
                            break;
                        }
                        coli++;
                    }
                    if(valid) cir.push_back(l);
                }
                l++;
            }


            //[> If we have less rows than columns, then transpose C_ij <]
            if(cic.size() <= cir.size()) {
                ci = cic;
            } else {
                ci = cir;
                C_ij = CT_ij;
                reversed = true;
            }

            if(ci.empty()) continue;

            int n_cij_before = C_ij.dim(1);


            // C_ij = CompCol_Mat_double(C_ij.dim(0), ci.size(), C_ij.NumNonzeros(),
            //         C_ij.val(MV_VecIndex()), C_ij.row_ind(MV_VecIndex()),
            //         C_ij.col_ptr(MV_VecIndex(0,ci.size())));
            // cout << C_ij << endl;
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
                    iv(k) = (double) -1.0;
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
    cout << "Size of C : " << nbcols - A.dim(1) << endl;
    size_c = nbcols - A.dim(1);
    n = nbcols;

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
