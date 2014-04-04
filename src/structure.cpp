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
    int guessPartitionsNumber = icntl[Controls::part_guess];

    if(guessPartitionsNumber == 1 && icntl[Controls::part_type] > 1){
        if (m_o == 1) {
            nbparts = 1;
        } else if (m_o <= 8) {
            nbparts = 2;
        } else if (m_o <= 1000) {
            nbparts = 4;
        } else if (m_o <= 50000) {
            nbparts = 8;
        } else if (m_o <= 100000) {
            nbparts = ceil((double)m_o / 10000);
        } else {
            nbparts = ceil((double)m_o / 20000);
        }
        cout << "Estimated number of partitions: " << nbparts  << endl;
        mpi::communicator world;
        parallel_cg =  nbparts < world.size() ? nbparts : world.size();
    }

    if (nbparts == 0){
        info[Controls::status] = -3;
        throw std::runtime_error("FATAL ERROR: Number of partitions is zero");
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
    if (nbparts == 1 && icntl[Controls::part_type] == 3) {
        cerr << "WARNING: PaToH is useless with a single partiton request, switching to automatic partitioning" << endl;
        icntl[Controls::part_type] = 2;
    }


    switch(icntl[Controls::part_type]){
        
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
            _imba   = dcntl[Controls::part_imbalance];
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
                info[Controls::status] = -4;
                stringstream err; 
                err <<"Error : PaToH Allocation problem : " << ret;
                throw std::runtime_error(err.str());
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
                    info[Controls::status] = -5;
                    throw std::runtime_error("FATAL ERROR: PaToH produced an empty partition, Try to reduce the imbalancing factor");
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
            info[Controls::status] = -6;
            throw std::runtime_error("Trying to use PaToH while it is not available");
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
        if(icntl[Controls::aug_type] == 0)
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
    if(icntl[Controls::aug_analysis] == 2){
        double f = 0;
        size_c = 1;
        cout << endl;
        while(size_c > 0 && f < 0.9){
            dcntl[Controls::aug_filter] = f;
            cout << "filter value : " << fixed << setprecision(5) << f << " gives : ";
            abcd::augmentMatrix(loc_parts);
            cout << endl << endl;
            f+=0.025;
        }
        exit(0);
    }

    if (icntl[Controls::aug_type] != 0) {
        abcd::augmentMatrix(loc_parts);
        cout << "   time to aug : " << MPI_Wtime() - t << endl;

        column_index.clear();
        for (unsigned int k = 0; k < (unsigned int)nbparts; k++) {

            CompCol_Mat_double part = loc_parts[k];

            // Build the column index of part
            std::vector<int> ci = getColumnIndex(
                                      part.colptr_ptr(), part.dim(1)
                                  );

            column_index.push_back(ci);
            ci_sizes.push_back(ci.size());

            parts[k] = CompRow_Mat_double(sub_matrix(part, ci));
        }
        if (icntl[Controls::aug_type] != 0)
            cout << "    time to part /w augmentation : " << MPI_Wtime() - t << endl;
        if (size_c == 0) {
            cerr << "WARNING: Size of C is zero, switching to classical cg" << endl;
            icntl[Controls::aug_type] = 0;
        }
    }
    // print only the size of C
    if(icntl[Controls::aug_analysis] == 1) exit(0);

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
    double filter_c = dcntl[Controls::aug_filter];
    stC = vector<int>(M.size(), -1);
    /*
     * Which augmentation to use:
     */
    if(icntl[Controls::aug_type] == 0){
        //[> No augmentation <]
        return;
    } else if (icntl[Controls::aug_type] == 1){
        /*
         * C_ij/-I augmentation
         */
        cijAugmentMatrix(M);
    } else if (icntl[Controls::aug_type] == 2){
        /*
         * A_ij/-A_ji augmentation
         */
        aijAugmentMatrix(M);

    } else {
        throw std::runtime_error("Unkown augmentation scheme.");
    }

}// [> -----  end of function abcd::augmentMatrix  ----- <]
