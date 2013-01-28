#include <abcd.h>
#include <iostream>
#include <fstream>

///@TODO Add the case of manual partitioning
void abcd::partitionMatrix()
{
    unsigned nbrows_per_part;
    unsigned row_sum = 0;

    if(nbparts == 0)
        throw - 13;
    if(nbparts < parallel_cg)
        throw - 101;

    nbrows_per_part = ceil(float(m)/float(nbparts));

    switch(partitioning_type){
        
        /*-----------------------------------------------------------------------------
         *  Uniform partitioning with a given nbrows
         *-----------------------------------------------------------------------------*/
        case 1:
            strow = VECTOR_int(nbparts);

            for(unsigned k = 0; k < nbparts; k++) {
                strow(k) = row_sum;
                row_sum += nbrows(k);
            }
            break;

        /*-----------------------------------------------------------------------------
         *  Uniform partitioning with only nbparts as input (generates nbrows)
         *-----------------------------------------------------------------------------*/
        case 2:
            strow = VECTOR_int(nbparts);
            nbrows = VECTOR_int(nbparts, nbrows_per_part);

            nbrows(nbparts - 1) = m - (nbparts - 1) * nbrows_per_part;

            for(unsigned k = 0; k < nbparts; k++) {
                strow(k) = row_sum;
                row_sum += nbrows(k);
            }
            break;
        /*-----------------------------------------------------------------------------
         *  PaToH partitioning
         *-----------------------------------------------------------------------------*/
        //case 3:

    }

}

void abcd::analyseFrame()
{

    std::vector<CompCol_Mat_double > loc_parts;
    std::vector<int> ci_sizes;


    double t  = MPI_Wtime();
    for(unsigned k = 0; k < nbparts; k++) {
        CompCol_Mat_double part = CSC_middleRows(A, strow[k], nbrows[k]);
        loc_parts.push_back(part);

        std::vector<int> ci;
        int j = 0;
        for(int i = 1; i <= loc_parts[k].dim(1); i++) {
            if(loc_parts[k].col_ptr(i) != loc_parts[k].col_ptr(i - 1))
                ci.push_back(j);
            j++;
        }
        column_index.push_back(ci);
    }
    //

    t= MPI_Wtime();
    abcd::augmentMatrix(loc_parts);
    cout << "time to aug : " << MPI_Wtime() -t << endl;
    exit(0);

    for(unsigned k = 0; k < nbparts; k++) {
        double t1, t2;
        // Build the column index of part
        column_index.clear();
        std::vector<int> ci;
        int j = 0;
        for(int i = 1; i <= loc_parts[k].dim(1); i++) {
            if(loc_parts[k].col_ptr(i) != loc_parts[k].col_ptr(i - 1))
                ci.push_back(j);
            j++;
        }
        column_index.push_back(ci);
        ci_sizes.push_back(ci.size());

        //int *last = std::unique(part.outerIndexPtr(), part.outerIndexPtr() + part.outerSize() + 1);
        //parts.push_back(SparseMatrix<double, RowMajor>(part.middleCols(0, ci.size())));
        VECTOR_int col_vect = loc_parts[k].t_col();
        //cout << col_vect << endl;
        int *last = std::unique(col_vect.t_vec(), col_vect.t_vec() + loc_parts[k].dim(1) + 1);

        partitions.push_back( 
                CompRow_Mat_double( 
                    CompCol_Mat_double( loc_parts[k].dim(0), ci.size(), loc_parts[k].NumNonzeros(),
                        loc_parts[k].t_val(), loc_parts[k].t_row(), col_vect(MV_VecIndex(0, ci.size())) 
                        ) 
                    )
                );
    }
    cout << "time to part : " << MPI_Wtime() -t << endl;

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::augmentMatrix
 *  Description:  Augments the matrix and build the C part in [A C]
 * =====================================================================================
 */
    void
abcd::augmentMatrix ( std::vector<CompCol_Mat_double> &M)
{
    double filter_c = dcntl[10];
    /*
     * Which augmentation to use:
     */
    if(icntl[10] == 0){
        //[> No augmentation <]
        return;
    } else if (icntl[10] == 1){
        /*
         * C_ij/-I augmentation
         */
        
        int nbcols = A.dim(1);
        int nz_c = 0;

        std::map<int,std::vector<CompCol_Mat_double> > C;
        std::map<int,std::vector<int> > stCols;

        //Eigen::SparseMatrix<double> Ti, Tj;
        for( int i = 0; i < M.size() - 1; i++ ){
            for ( int j = i+1; j < M.size(); j++ ) {
                std::vector<int> intersect;
                std::set_intersection(column_index[i].begin(), column_index[i].end(),
                                        column_index[j].begin(), column_index[j].end(),
                                        std::back_inserter(intersect));
                if (intersect.size() == 0) continue;

                CompCol_Mat_double C_ij;
                {
                    CompRow_Mat_double A_ij = CompRow_Mat_double(sub_matrix(M[i], intersect));
                    CompRow_Mat_double A_ji = CompRow_Mat_double(sub_matrix(M[j], intersect));
                    C_ij = smmtm(A_ij, A_ji);
                }

                if(C_ij.NumNonzeros() == 0) continue;
                //cout << C_ij.NumNonzeros() << endl;

                /*-----------------------------------------------------------------------------
                 *  Compress C_ij
                 *-----------------------------------------------------------------------------*/
                std::vector<int> cic, cir, ci;
                bool reversed = false;

                //[> Compute the Column compression <]
                int l = 0;
                for(int k = 1; k <= C_ij.dim(1); k++) {
                    if(C_ij.col_ptr(k) != C_ij.col_ptr(k - 1)){
                        bool valid = false;
                        int coli = C_ij.col_ptr(k-1);
                        while(coli < C_ij.col_ptr(k)){
                            if(abs(C_ij.val(coli)) >= filter_c){
                                valid = true;
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
                            if(abs(CT_ij.val(coli)) >= filter_c){
                                valid = true;
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

                if(ci.size() == 0) continue;

                int n_cij_before = C_ij.dim(1);

                int *last = std::unique(C_ij.colptr_.t_vec(),
                        C_ij.colptr_.t_vec() + C_ij.dim(1) + 1);


                C_ij = CompCol_Mat_double(C_ij.dim(0), ci.size(), C_ij.NumNonzeros(),
                        C_ij.val(MV_VecIndex()), C_ij.row_ind(MV_VecIndex()),
                        C_ij.col_ptr(MV_VecIndex(0,ci.size())));


                //// Build compressed I
                CompCol_Mat_double I;
                {
                    VECTOR_int ir(C_ij.dim(1));
                    VECTOR_int ic(C_ij.dim(1)+1);
                    VECTOR_double iv(C_ij.dim(1), -1);
                    for(int k = 0; k < C_ij.dim(1); k++){
                        ir(k) = ci[k];
                        //ir(k) = k;
                        ic(k) = k;
                    }
                    ic(C_ij.dim(1)) = C_ij.dim(1);

                    I = CompCol_Mat_double(n_cij_before, C_ij.dim(1), C_ij.dim(1),
                            iv, ir, ic);
                }


                if(!reversed){
                    stCols[i].push_back(nbcols);
                    stCols[j].push_back(nbcols);
                    C[i].push_back(C_ij);
                    C[j].push_back(I);
                } else {
                    stCols[i].push_back(nbcols);
                    stCols[j].push_back(nbcols);
                    C[i].push_back(I);
                    C[j].push_back(C_ij);
                }

                nz_c += C_ij.NumNonzeros() + I.NumNonzeros();
                nbcols += C_ij.dim(1);
            }

        }
        cout << "Size of C : " << nbcols - A.dim(1) << endl;

        // Augment the matrices
        for(int k = 0; k < M.size(); k++){
            // now augment each partition!
            M[k] = concat_columns(M[k], C[k], stCols[k]);
            M[k] = resize_columns(M[k], nbcols);
        }


    } else if (icntl[10] == 2){
        /*
         * A_ij/-A_ji augmentation
         */
        int nbcols = A.dim(1);
        std::map<int,std::vector<CompCol_Mat_double> > C;
        std::map<int,std::vector<int> > stCols;


        for( int i = 0; i < M.size() - 1; i++ ){
            for ( int j = i+1; j < M.size(); j++ ) {
                std::vector<int> intersect;
                std::set_intersection(column_index[i].begin(), column_index[i].end(),
                                        column_index[j].begin(), column_index[j].end(),
                                        std::back_inserter(intersect));
                if (intersect.size() == 0) continue;

                CompCol_Mat_double A_ij = sub_matrix(M[i], intersect);
                CompCol_Mat_double A_ji = sub_matrix(M[j], intersect);
                for(int k = 0; k < A_ji.NumNonzeros(); k++)
                    A_ji.val(k) *= double(-1);


                if(filter_c != 0) {
                    std::vector<int> selected_cols;
                    std::vector<double> frob_ij, mu;
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

                    nu = (frob_sum / A_ij.dim(1)) / sqrt(card_max);

                    for (int k = 0; k < A_ij.dim(1); ++k){
                        VECTOR_int A_ij_k_ind, A_ji_k_ind;
                        VECTOR_double A_ij_k = middleCol(A_ij, k, A_ij_k_ind);
                        VECTOR_double A_ji_k = middleCol(A_ji, k, A_ji_k_ind);

                        double inf_ij = infNorm(A_ij_k);
                        double inf_ji = infNorm(A_ji_k);

                        double p = 0, q = 0;
                        for (int l = 0; l < A_ij_k_ind.size(); ++l){
                            if (abs(A_ij_k(l)) >= nu/inf_ji) p++;
                        }
                        for (int l = 0; l < A_ji_k_ind.size(); ++l){
                            if (abs(A_ji_k(l)) >= nu/inf_ij) q++;
                        }

                        p = ( p==0 ? A_ij_k_ind.size() : p );
                        q = ( q==0 ? A_ji_k_ind.size() : q );

                        double mu_ij_k = frob_ij[k] / sqrt(p*q);
                        mu.push_back(mu_ij_k);

                        if(mu_ij_k >= filter_c){
                            selected_cols.push_back(k);
                        }
                    }

                    if (selected_cols.size() == 0) continue;

                    A_ij = sub_matrix(A_ij, selected_cols);
                    A_ji = sub_matrix(A_ji, selected_cols);
                }


                stCols[i].push_back(nbcols);
                stCols[j].push_back(nbcols);
                C[i].push_back(A_ij);
                C[j].push_back(A_ji);

                nbcols += A_ij.dim(1);
            }
        }

        cout << "Size of C : " << nbcols - A.dim(1) << endl;


        // Augment the matrices
        for(int k = 0; k < M.size(); k++){
            // now augment each partition!
            M[k] = concat_columns(M[k], C[k], stCols[k]);
            M[k] = resize_columns(M[k], nbcols);
        }

    } else if (icntl[10] == 3){
        /*
         * SVD augmentation
         */

    }

}// [> -----  end of function abcd::augmentMatrix  ----- <]
