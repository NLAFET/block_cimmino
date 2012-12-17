#include <abcd.h>

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

            nbrows(nbparts - 1) += m;
            for (int i = 0; i < nbrows.size(); i++)
                nbrows(nbparts - 1) -= nbrows(i);

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
    int *ip, *jp;
    double *vp;

    int st_p = 0, ed_pt;
    int sssm = 0;
    CompRow_Mat_double part_r;
    double t  = MPI_Wtime();

    VECTOR_int ar = A.t_row();
    VECTOR_int ac = A.t_col();
    VECTOR_double av = A.t_val();
    int ind_st, ind_ed, rsp;
    for(unsigned k = 0; k < nbparts; k++) {
        ind_st = A.row_ptr(strow[k]);
        ind_ed = A.row_ptr(strow[k] + nbrows[k]) - 1;
        rsp = ar(strow[k]);

        VECTOR_int rpt(nbrows[k]+1);
        rpt = ar(MV_VecIndex(strow[k], strow[k]+nbrows[k]));
        for(int i = 0; i <= nbrows[k]; i++){
            rpt[i] -= rsp;
        }

        // Our k-th partition
        part_r = CompRow_Mat_double (
                nbrows[k], n, ind_ed - ind_st,
                av(MV_VecIndex(ind_st, ind_ed)),
                rpt,
                ac(MV_VecIndex(ind_st, ind_ed))
                );

        loc_parts.push_back(CompCol_Mat_double(part_r));
    }
    cout << "time to part : " << MPI_Wtime() -t << endl;

    double t= MPI_Wtime();
    abcd::augmentMatrix(loc_parts);
    cout << "time to aug : " << MPI_Wtime() -t << endl;

    //for(unsigned k = 0; k < nbparts; k++) {
        //double t1, t2;
        //// Build the column index of part
        //std::vector<int> ci;
        //int j = 0;
        //for(int i = 1; i <= loc_parts[k].outerSize(); i++) {
            //if(loc_parts[k].outerIndexPtr()[i] != loc_parts[k].outerIndexPtr()[i - 1])
                //ci.push_back(j);
            //j++;
        //}
        //column_index.push_back(ci);
        //ci_sizes.push_back(ci.size());

        ////int *last = std::unique(part.outerIndexPtr(), part.outerIndexPtr() + part.outerSize() + 1);
        ////parts.push_back(SparseMatrix<double, RowMajor>(part.middleCols(0, ci.size())));
        //int *last = std::unique(loc_parts[k].outerIndexPtr(),
                //loc_parts[k].outerIndexPtr() + loc_parts[k].outerSize() + 1);
        //parts.push_back(SparseMatrix<double, RowMajor>(loc_parts[k].middleCols(0, ci.size())));
        ////loc_parts[k] = NULL;
    //}

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::augmentMatrix
 *  Description:  Augments the matrix and build the C part in [A C]
 * =====================================================================================
 */
    void
abcd::augmentMatrix ( std::vector<CompCol_Mat_double> > &M)
{
    /*
     * Which augmentation to use:
     */
    if(icntl[10] == 0){
        /* No augmentation */
        return;
    } else if (icntl[10] == 1){
        /*
         * C_ij/-I augmentation
         */
        
        int nbcols = mtx.cols();
        Eigen::SparseMatrix<double> Ti, Tj;
        for( int i = 0; i < M.size() - 1; i++ ){
            for ( int j = i+1; j < M.size(); j++ ) {
                // resize A_ij and A_ji as Ti and Tj respectively
                Ti.resize(M[i].rows(), nbcols + M[j].rows());
                Tj.resize(M[j].rows(), nbcols + M[j].rows());

                // put the old A_ij and A_ji at the begining
                Ti.middleCols(0, M[i].cols()) = M[i];
                Tj.middleCols(0, M[j].cols()) = M[j];

                // Compute C_ij, The added part should not disturb the
                // computation
                Eigen::SparseMatrix<double> C_ij = Ti*Tj.transpose();
                if(C_ij.nonZeros() == 0) continue;

                /*-----------------------------------------------------------------------------
                 *  Compress C_ij
                 *-----------------------------------------------------------------------------*/
                std::vector<int> cic, cir, ci;
                bool reversed = false;

                /* Compute the Column compression */
                int l = 0;
                for(int k = 1; k <= C_ij.outerSize(); k++) {
                    if(C_ij.outerIndexPtr()[k] != C_ij.outerIndexPtr()[k - 1])
                        cic.push_back(l);
                    l++;
                }
                /* Compute the Row compression */
                Eigen::SparseMatrix<double> CT_ij = C_ij.transpose();
                l = 0;
                for(int k = 1; k <= CT_ij.outerSize(); k++) {
                    if(CT_ij.outerIndexPtr()[k] != CT_ij.outerIndexPtr()[k - 1])
                        cir.push_back(l);
                    l++;
                }

                /* If we have less rows than columns, then transpose C_ij */
                if(cic.size() <= cir.size()) {
                    ci = cic;
                } else {
                    ci = cir;
                    C_ij = CT_ij;
                    reversed = true;
                }

                int *last = std::unique(C_ij.outerIndexPtr(),
                        C_ij.outerIndexPtr() + C_ij.outerSize() + 1);
                C_ij = Eigen::SparseMatrix<double>(C_ij.middleCols(0, ci.size()));

                // Build compressed I
                Eigen::SparseMatrix<double> I(M[i].rows(), C_ij.cols());
                for(int k = 0; k < C_ij.cols(); k++)
                    I.coeffRef(ci[k],k) = -1;


                if(!reversed){
                    // Add C_ij to A_ij
                    Ti.middleCols(nbcols, C_ij.cols()) = C_ij;
                    // Add I to A_ji
                    Tj.middleCols(nbcols, C_ij.cols()) = I;

                    //Compress Ti and Tj as we over-estimated their size
                    Ti = Eigen::SparseMatrix<double>(Ti.middleCols(0, nbcols + C_ij.cols()));
                    Tj = Eigen::SparseMatrix<double>(Tj.middleCols(0, nbcols + C_ij.cols()));
                } else {
                    // Add I to A_ij
                    Ti.middleCols(nbcols, C_ij.cols()) = I;
                    // Add C_ij^T to A_ji
                    Tj.middleCols(nbcols, C_ij.cols()) = C_ij;

                    //Compress Ti and Tj as we over-estimated their size
                    Ti = Eigen::SparseMatrix<double>(Ti.middleCols(0, nbcols + C_ij.cols()));
                    Tj = Eigen::SparseMatrix<double>(Tj.middleCols(0, nbcols + C_ij.cols()));
                }

                M[i] = Ti;
                M[j] = Tj;


                nbcols += C_ij.cols();
            }

        }
        cout << "Size of C : " << nbcols - mtx.cols() << endl;
        //exit(0);

    } else if (icntl[10] == 2){
        /*
         * A_ij/-A_ji augmentation
         */

    } else if (icntl[10] == 3){
        /*
         * SVD augmentation
         */

    }

}		/* -----  end of function abcd::augmentMatrix  ----- */
