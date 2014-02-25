#include <abcd.h>
#include "mat_utils.h"
#include <iostream>
#include <fstream>
#include <patoh.h>
#include <vect_utils.h>
#include <boost/lambda/lambda.hpp>

using namespace boost::lambda;

void abcd::partitionMatrix()
{
    unsigned handled_rows = 0;
    unsigned ceil_per_part, floor_per_part;
    unsigned row_sum = 0;

    if (nbparts == 0){
        cerr << "FATAL ERROR: Number of partitions is zero" << endl;
        throw -1;
    }
    if (nbparts < parallel_cg) {
        cerr << "ERROR: Number of partitions is smaller than the number of parallel_cg" << endl;
        cerr << "WARNING: Increasing the number of partitions from " << nbparts
            << " up to " << parallel_cg << endl;
        nbparts = parallel_cg;
    }
    if(nbparts > m) {
        cerr << "ERROR: Number of partitions is larger than the number of rows" << endl;
        cerr << "WARNING: Decreasing the number of partitions from " << nbparts
            << " down to " << m << endl;
        nbparts = parallel_cg;
    }
    if (nbparts == 1 && partitioning_type == 3) {
        cerr << "WARNING: PaToH is useless with a single partiton request, switching to automatic partitioning" << endl;
        partitioning_type = 2;
    }

    if(guessPartitionsNumber == 1 && partitioning_type > 1){
        // if the number of rows is less than 1k
        if (m_o <= 1000) {
            nbparts = 4;
        }
        // if the number of rows is less than 10k
        else if (m_o <= 10000) {
            nbparts = 8;
        // if the number of rows is between 10k and 160k
        } else if (m_o <= 160000) {
            nbparts = 16;
        // if the number of rows is larger than 100k
        } else if (m_o <= 1000000) {
            nbparts = ceil(m_o / 15000);
        }
        cout << "Estimated number of partitions: " << nbparts  << endl;
    }

    switch(partitioning_type){
        
        /*-----------------------------------------------------------------------------
         *  Uniform partitioning with a given nbrows
         *-----------------------------------------------------------------------------*/
        case 1:
            strow = VECTOR_int(nbparts);

            for(unsigned int k = 0; k < (unsigned int)nbparts; k++) {
                strow(k) = row_sum;
                row_sum += nbrows(k);
            }
            break;

        /*-----------------------------------------------------------------------------
         *  Uniform partitioning with only nbparts as input (generates nbrows)
         *-----------------------------------------------------------------------------*/
        case 2:
            ceil_per_part = ceil(float(m_o)/float(nbparts));
            floor_per_part = floor(float(m_o)/float(nbparts));

            strow = VECTOR_int(nbparts);
            nbrows = VECTOR_int(nbparts);

            // alternate the number of rows, they will not be that equal
            // but at least we will have simmilar number of row
            for(unsigned k = 0; k < (unsigned) nbparts; k+=2) {
                nbrows(k) = ceil_per_part;
                handled_rows += ceil_per_part;
            }
            for(unsigned k = 1; k < (unsigned) nbparts; k+=2) {
                nbrows(k) = floor_per_part;
                handled_rows += floor_per_part;
            }

            nbrows(nbparts - 1) += m_o - handled_rows;

            for(unsigned k = 0; k < (unsigned)nbparts; k++) {
                strow(k) = row_sum;
                row_sum += nbrows(k);
            }
            break;
        /*-----------------------------------------------------------------------------
         *  PaToH partitioning
         *-----------------------------------------------------------------------------*/
        case 3:
#ifdef PATOH
            PaToH_Parameters args;
            int _c, _n, _nconst, _imba, _ne, *cwghts, *nwghts, *xpins, *pins, *partvec,
                cut, *partweights, ret;
            char cutdef[] = "CUT";

            CompCol_Mat_double t_A = Coord_Mat_double(A);

            double t = MPI_Wtime();
            cout << "[-] launching PaToH" << endl;

            PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_DEFAULT);
            args._k = nbparts;
            _c = m_o;
            _n = n_o;
            _nconst = 1;
            _imba   = dcntl[8];
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
                throw -1;
            }


            args.final_imbal    = _imba;
            args.init_imbal     = _imba * 2.0;
            args.seed           = 1;

            //args.initp_alg      = 12;

            partvec     = new int[_c];
            partweights = new int[args._k * _nconst];

            PaToH_Part(&args, _c, _n, _nconst, 0, cwghts, nwghts,
                            xpins, pins, NULL, partvec, partweights, &cut);
            for (int i = 0; i < nbparts; i++) {
                if (partweights[i] == 0) {
                    cerr << "FATAL ERROR: PaToH produced an empty partition" << endl
                        << "Try to reduce the imbalancing factor" << endl;
                    throw -1;
                }
            }

            row_perm = sort_indexes(partvec, _c);

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
            cout << "    Done with PaToH, time : " << MPI_Wtime() - t << endl;
            t = MPI_Wtime();

            A = CompRow_Mat_double(m_o, n_o, nz_o, val, ir, jc);

            nbrows = VECTOR_int(partweights, nbparts);
            strow = VECTOR_int(nbparts);

            for(unsigned k = 0; k < nbparts; k++) {
                strow(k) = row_sum;
                row_sum += nbrows(k);
            }

           if(write_problem.length() != 0) {
                ofstream f;
                f.open(write_problem.c_str());
                f << "%%MatrixMarket matrix coordinate real general\n";
                f << A.dim(0) << " " << A.dim(1) << " " << A.NumNonzeros() << "\n";
                for(int i = 0; i < m_o; i++){
                    for(int j = ir[i]; j< ir[i + 1]; j++){
                        f << i + 1 << " " << jc[j] + 1 << " " << val[j] << "\n";
                    }
                }
                f.close();

            }


            cout << "    Finished Partitioning, time : " << MPI_Wtime() - t << endl;
            delete[] ir, jc, val, partvec, partweights, cwghts, pins, xpins, nwghts,
                ir, jc, val;
            PaToH_Free();
            t_A = CompCol_Mat_double();
#else
	    throw -99;
#endif
            break;
    }

    if(write_problem.length() != 0) {
        string parts = write_problem + "_parts";
        ofstream f;
        f.open(parts.c_str());

        for(unsigned int k = 0; k < (unsigned int)nbparts; k++) {
            f << nbrows[k] << "\n";
        }

        f.close();
    }
}

void abcd::analyseFrame()
{

    std::vector<CompCol_Mat_double > loc_parts;
    loc_parts.reserve(nbparts);
    std::vector<int> ci_sizes;


    double t  = MPI_Wtime();

    column_index.reserve(nbparts);
    cout << "[+] Creating partitions"<< flush;
    
    for (unsigned int k = 0; k < (unsigned int)nbparts; k++) {
        CompCol_Mat_double part = CSC_middleRows(A, strow[k], nbrows[k]);
        
        int *col_ptr = part.colptr_ptr();
        std::vector<int> ci = getColumnIndex(col_ptr, part.dim(1));
        column_index.push_back( ci );

        // if no augmentation, then create the parts
        if(icntl[10] == 0)
        {
            parts[k] = CompRow_Mat_double(sub_matrix(part, ci));
        } else 
        {
            loc_parts.push_back(part);
        }
    }
    cout << ", done in " << MPI_Wtime() - t<< endl;
    //
    t= MPI_Wtime();

    // test augmentation!
    if(icntl[11] == 2){
        double f = 0;
        size_c = 1;
        cout << endl;
        while(size_c > 0 && f < 0.9){
            dcntl[10] = f;
            cout << "filter value : " << fixed << setprecision(5) << f << " gives : ";
            abcd::augmentMatrix(loc_parts);
            cout << endl << endl;
            f+=0.025;
        }
        exit(0);
    }

    if (icntl[10] != 0) {
        abcd::augmentMatrix(loc_parts);
        cout << "   time to aug : " << MPI_Wtime() - t << endl;

        column_index.clear();
        for (unsigned int k = 0; k < (unsigned int)nbparts; k++) {

            CompCol_Mat_double part = loc_parts[k];
            //int *col_vect_ptr = part.colptr_ptr();

            // Build the column index of part
            std::vector<int> ci = getColumnIndex(
                                      part.colptr_ptr(), part.dim(1)
                                  );

            column_index.push_back(ci);
            ci_sizes.push_back(ci.size());

            //parts[k] =
                //CompRow_Mat_double(
                    //CompCol_Mat_double(part.dim(0), ci.size(),
                                       //part.NumNonzeros(),
                                       //part.val_ptr(),
                                       //part.rowind_ptr(),
                                       //col_vect_ptr
                                      //)
                //);
            parts[k] = CompRow_Mat_double(sub_matrix(part, ci));
        }
        if (icntl[10] != 0)
            cout << "    time to part /w augmentation : " << MPI_Wtime() - t << endl;
        if (size_c == 0) {
            cerr << "WARNING: Size of C is zero, switching to classical cg" << endl;
            icntl[10] = 0;
        }
    }
    // print only the size of C
    if(icntl[11] == 1) exit(0);

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
    stC = vector<int>(M.size(), -1);
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
        for( size_t i = 0; i < M.size() - 1; i++ ){
            for ( size_t j = i+1; j < M.size(); j++ ) {
                std::vector<int> intersect;
                std::set_intersection(column_index[i].begin(), column_index[i].end(),
                                        column_index[j].begin(), column_index[j].end(),
                                        std::back_inserter(intersect));
                if (intersect.size() == 0) continue;

                CompCol_Mat_double C_ij;
                {
                    CompRow_Mat_double A_ij = CompRow_Mat_double(sub_matrix(M[i], intersect));
                    CompRow_Mat_double A_ji = CompRow_Mat_double(sub_matrix(M[j], intersect));
                    CompRow_Mat_double A_jiT = csr_transpose(A_ji);
                    C_ij = spmm(A_ij, A_jiT);
                }

                if(C_ij.NumNonzeros() == 0) continue;
                //cout << C_ij.NumNonzeros() << endl;

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
                            if(icntl[15] != 0){ 
                                if(abs(C_ij.val(coli)) >= dcntl[15]){
                                    selected_S_columns.push_back( nbcols + k - n_o);
                                } else {
                                    skipped_S_columns.push_back( nbcols + k - n_o);
                                }
                            }

                            if(abs(C_ij.val(coli)) >= filter_c){
                                valid = true;
                                break;
                            }

                            if( icntl[15] != 2 ){ // don't reduce, we just need the selected columns!
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
                            if(icntl[15] != 0){ 
                                if(abs(CT_ij.val(coli)) >= dcntl[15]){
                                    selected_S_columns.push_back( nbcols + k - n_o);
                                } else {
                                    skipped_S_columns.push_back( nbcols + k - n_o);
                                }
                            }

                            if(abs(CT_ij.val(coli)) >= filter_c){
                                valid = true;
                                break;
                            }

                            if( icntl[15] != 2 ){ // don't reduce, we just need the selected columns!
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

                if(ci.size() == 0) continue;

                int n_cij_before = C_ij.dim(1);


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

        if(icntl[11] != 0) return;

        // Augment the matrices
        for(size_t k = 0; k < M.size(); k++){
            if(stCols[k].size() == 0) continue;
            // now augment each partition!
            stC.push_back(stCols[k][0]);
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


        for( size_t i = 0; i < M.size() - 1; i++ ){
            for ( size_t j = i+1; j < M.size(); j++ ) {
                std::vector<int> intersect;
                std::set_intersection(column_index[i].begin(), column_index[i].end(),
                                        column_index[j].begin(), column_index[j].end(),
                                        std::back_inserter(intersect));
                if (intersect.size() == 0) continue;

                CompCol_Mat_double A_ij = sub_matrix(M[i], intersect);
                CompCol_Mat_double A_ji = sub_matrix(M[j], intersect);

                for(int k = 0; k < A_ji.NumNonzeros(); k++)
                    A_ji.val(k) *= double(-1);


                if(filter_c != 0 || dcntl[15] != 0) {
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

                        if(icntl[15] != 0){ 
                            if (dcntl[15] < 0) {
                                if ((nbcols + k - n_o) % abs((int)dcntl[15]) == 0)
                                    selected_S_columns.push_back( nbcols + k - n_o);
                                else
                                    skipped_S_columns.push_back( nbcols + k - n_o);
                            } else {
                                if (mu_ij_k >= dcntl[15])
                                    selected_S_columns.push_back( nbcols + k - n_o);
                                else
                                    skipped_S_columns.push_back( nbcols + k - n_o);
                            }
                        }

                    }

                    if (selected_cols.size() == 0) continue;

                    if( icntl[15] != 2 ) { // don't reduce the A_ij/A_ji, we just need the selected columns!
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

        cout << "Size of C : " << nbcols - A.dim(1) << endl;
        size_c = nbcols - A.dim(1);
        n = nbcols;

        if(icntl[11] != 0) return;

        // Augment the matrices
        for(size_t k = 0; k < M.size(); k++){
            if(stCols[k].size() == 0) continue;
            // now augment each partition!
            stC[k] = stCols[k][0];
            M[k] = concat_columns(M[k], C[k], stCols[k]);
            M[k] = resize_columns(M[k], nbcols);
        }

    } else if (icntl[10] == 3){
        /*
         * SVD augmentation
         */

    }

}// [> -----  end of function abcd::augmentMatrix  ----- <]
