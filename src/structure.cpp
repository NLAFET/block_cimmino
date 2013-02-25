#include <abcd.h>
#include <iostream>
#include <fstream>
#include <patoh.h>
#include <boost/lambda/lambda.hpp>

using namespace boost::lambda;

///@TODO Add the case of manual partitioning
void abcd::partitionMatrix()
{
    unsigned nbrows_per_part;
    unsigned row_sum = 0;

    if(nbparts == 0)
        throw - 13;
    if(nbparts < parallel_cg)
        throw - 101;

    nbrows_per_part = ceil(float(m_o)/float(nbparts));

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

            nbrows(nbparts - 1) = m_o - (nbparts - 1) * nbrows_per_part;

            for(unsigned k = 0; k < nbparts; k++) {
                strow(k) = row_sum;
                row_sum += nbrows(k);
            }
            break;
        /*-----------------------------------------------------------------------------
         *  PaToH partitioning
         *-----------------------------------------------------------------------------*/
        case 3:
            PaToH_Parameters args;
            int _c, _n, _nconst, _imba, _ne, *cwghts, *nwghts, *xpins, *pins, *partvec,
                cut, *partweights, ret;
            char cutdef[] = "CUT";

            CompCol_Mat_double t_A = A;

            PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_DEFAULT);
            args._k = nbparts;
            _c = m_o;
            _n = n_o;
            _nconst = 1;
            _imba   = 0.03;
            _ne     = nz_o;

            //xpins   = t_A.colptr_ptr();
            //pins    = t_A.rowind_ptr();
            xpins   = new int[_n + 1];
            pins    = new int[nz_o];

            for(int i = 0; i <= _n; i++){
                xpins[i] = t_A.col_ptr(i);
            }
            for(int i = 0; i < _ne; i++){
                pins[i] = t_A.row_ind(i);
            }

            cwghts  = new int[_c*_nconst];
            //using boost lambdas
            for_each(cwghts, cwghts + _c, _1 = 1);

            nwghts  = NULL;

            if( ret = PaToH_Alloc(&args, _c, _n, _nconst, cwghts, nwghts, xpins, pins) ){
                cerr << "Error : PaToH Allocation problem : " << ret << endl;
                throw -102;
            }

            args.final_imbal    = _imba;
            args.init_imbal     = _imba * 2.0;
            args.seed           = 0;

            partvec     = new int[_c];
            partweights = new int[args._k * _nconst];

            PaToH_Part(&args, _c, _n, _nconst, 0, cwghts, nwghts,
                            xpins, pins, NULL, partvec, partweights, &cut);

            std::vector<int> row_perm = sort_indexes(partvec, _c);

            // Permutation
            int *iro = A.rowptr_ptr();
            int *jco = A.colind_ptr();
            double *valo = A.val_ptr();

            int *ir = new int[m_o + 1];
            int *jc = new int[nz_o];
            double *val = new double[nz_o];

            int sr = 0;
            for(int i = 0; i < m_o; i++){
                int cur = row_perm[i];
                ir[i] = sr;
                for(int j = 0; j < iro[cur+1] - iro[cur]; j++){
                    jc[ir[i] + j] = jco[iro[cur] + j];
                    val[ir[i] + j] = valo[iro[cur] + j];
                }
                sr += iro[cur + 1] - iro[cur];
            }
            ir[m_o] = nz_o;


            A = CompRow_Mat_double(m_o, n_o, nz_o, val, ir, jc);

            if(write_problem.length() != 0) {
                ofstream f;
                f.open(write_problem.c_str());
                f << "%%MatrixMarket matrix coordinate real general\n";
                f << A.dim(0) << " " << A.dim(1) << " " << A.NumNonzeros() << "\n";
                for(int i = 0; i < m_o; i++){
                    for(int j = ir[i]; j< ir[i + 1]; j++){
                        f << row_perm[i] << " " << i + 1 << " " << jc[j] + 1 << " " << val[j] << "\n";
                    }
                }
                f.close();
            }

            nbrows = VECTOR_int(partweights, nbparts);
            strow = VECTOR_int(nbparts);

            for(unsigned k = 0; k < nbparts; k++) {
                strow(k) = row_sum;
                row_sum += nbrows(k);
            }

            cout << "Done partitioning" << endl;
            break;
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

    if(icntl[11] == 1) exit(0);
    if(icntl[11] == 2){
        double f = 0;
        size_c = 1;
        cout << endl;
        while(size_c > 0 && f < 0.9){
            dcntl[10] = f;
            cout << "filter value : " << fixed << setprecision(3) << f << " gives : ";
            abcd::augmentMatrix(loc_parts);
            cout << endl << endl;
            f+=25.0/1000;
        }
        exit(0);
    }

    abcd::augmentMatrix(loc_parts);
    cout << "time to aug : " << MPI_Wtime() -t << endl;

    column_index.clear();
    for(unsigned k = 0; k < nbparts; k++) {
        double t1, t2;
        // Build the column index of part
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
        int *last = std::unique(col_vect.ptr(), col_vect.ptr() + loc_parts[k].dim(1) + 1);

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
                    C_ij = spmm(A_ij, A_ji);
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

                int *last = std::unique(C_ij.colptr_ptr(),
                        C_ij.colptr_ptr() + C_ij.dim(1) + 1);


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
        size_c = nbcols - A.dim(1);
        n = nbcols;

        if(icntl[11] != 1) return;

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
        size_c = nbcols - A.dim(1);
        n = nbcols;

        if(icntl[11] != 1) return;

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
