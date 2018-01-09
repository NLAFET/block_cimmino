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

#include <abcd.h>
#include "mat_utils.h"
#include <iostream>
#include <fstream>
#include <vect_utils.h>
#include <boost/lambda/lambda.hpp>

#ifdef PATOH
#include <patoh.h>
#endif //PATOH

using namespace boost::lambda;

void abcd::partitionMatrix()
{
    unsigned handled_rows = 0;
    unsigned ceil_per_part, floor_per_part;
    unsigned row_sum = 0;
    int guessPartitionsNumber = icntl[Controls::part_guess];

    if(guessPartitionsNumber == 1 && icntl[Controls::part_type] > 1){
        if (m_o == 1) {
            icntl[Controls::nbparts] = 1;
        } else if (m_o <= 8) {
            icntl[Controls::nbparts] = 2;
        } else if (m_o <= 1000) {
            icntl[Controls::nbparts] = 4;
        } else if (m_o <= 50000) {
            icntl[Controls::nbparts] = 8;
        } else if (m_o <= 100000) {
            icntl[Controls::nbparts] = ceil((double)m_o / 10000);
        } else {
            icntl[Controls::nbparts] = ceil((double)m_o / 20000);
        }
        LINFO << "Estimated number of partitions: " << icntl[Controls::nbparts];
        parallel_cg =  icntl[Controls::nbparts] < comm.size() ? icntl[Controls::nbparts] : comm.size();
    }

    if (icntl[Controls::nbparts] == 0){
        info[Controls::status] = -3;
        mpi::broadcast(comm, info[Controls::status], 0);
        throw std::runtime_error("FATAL ERROR: Number of partitions is zero");
    }
    if (icntl[Controls::nbparts] < parallel_cg) {
        // the user is not supposed to know about the parallel_cg
        // LERROR << "ERROR: Number of partitions is smaller than the number of parallel_cg";
        LERROR << "Oops! This should not happen";
        
        LWARNING << "WARNING: Increasing the number of partitions from " << icntl[Controls::nbparts]
                 << " up to " << parallel_cg;
        
        icntl[Controls::nbparts] = parallel_cg;
    }
    if(icntl[Controls::nbparts] > m) {
        LERROR << "ERROR: Number of partitions is larger than the number of rows";
        
        LWARNING << "WARNING: Decreasing the number of partitions from " << icntl[Controls::nbparts]
                 << " down to " << m;
        
        icntl[Controls::nbparts] = m;
    }

    if (icntl[Controls::nbparts] == 1 && icntl[Controls::part_type] == 3) {
        LWARNING << "WARNING: PaToH is useless with a single partiton request, switching to automatic partitioning";
        icntl[Controls::part_type] = 2;
    }

    if (icntl[Controls::part_type] == 1 && nbrows.size() == 0) {
        info[Controls::status] = -4;
        mpi::broadcast(comm, info[Controls::status], 0);
        throw std::runtime_error("nbrows not initialized for manual partitioning");
    }
    


    switch(icntl[Controls::part_type]){
        
        /*-----------------------------------------------------------------------------
         *  Uniform partitioning with a given nbrows
         *-----------------------------------------------------------------------------*/
    case 1:
        strow.resize(icntl[Controls::nbparts]);

        for(unsigned int k = 0; k < (unsigned int)icntl[Controls::nbparts]; k++) {
            strow[k] = row_sum;
            row_sum += nbrows[k];
            if (nbrows[k] == 0) {
                info[Controls::status] = -6;
                mpi::broadcast(comm, info[Controls::status], 0);
                throw std::runtime_error("FATAL ERROR: You requested an empty partition.");
            }
        }
        break;

        /*-----------------------------------------------------------------------------
         *  Uniform partitioning with only icntl[Controls::nbparts] as input (generates nbrows)
         *-----------------------------------------------------------------------------*/
    case 2:
	{
	        unsigned floor_per_part;
	        floor_per_part = floor(float(m_o)/float(icntl[Controls::nbparts]));
	        int remain = m_o - ( floor_per_part * icntl[Controls::nbparts]);
	        strow =  std::vector<int>(icntl[Controls::nbparts]);
	        nbrows = std::vector<int>(icntl[Controls::nbparts]);
	 
	        for(unsigned k = 0; k < (unsigned) icntl[Controls::nbparts]; k++) {
	           int cnt = floor_per_part;
	           if(remain >0){
	               remain--;
	               cnt++;
	           }
	           nbrows[k] = cnt;
	        }
	 
	        for(unsigned k = 0; k < (unsigned)icntl[Controls::nbparts]; k++) {
	            strow[k] = row_sum;
	            row_sum += nbrows[k];
	        }
	}
        break;
        /*-----------------------------------------------------------------------------
         *  PaToH partitioning
         *-----------------------------------------------------------------------------*/
    case 3:
	{
#ifdef PATOH
        PaToH_Parameters args;
        int _c, _n, _nconst, _ne, *cwghts, *nwghts, *xpins, *pins, *partvec,
            cut, *partweights, ret;
	double _imba;

        CompCol_Mat_double t_A = Coord_Mat_double(A);

        double t = MPI_Wtime();
        LINFO << "Launching PaToH";

        PaToH_Initialize_Parameters(&args, PATOH_CUTPART, PATOH_SUGPARAM_DEFAULT);
        args._k = icntl[Controls::nbparts];
        _c = m_o;
        _n = n_o;
        _nconst = 1;
        _imba   = dcntl[Controls::part_imbalance];
        _ne     = nz_o;

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
            info[Controls::status] = -5;
            mpi::broadcast(comm, info[Controls::status], 0);
            stringstream err; 
            LERROR <<"Error : PaToH Allocation problem : " << ret;
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
        for (int i = 0; i < icntl[Controls::nbparts]; i++) {
            if (partweights[i] == 0) {
                info[Controls::status] = -6;
                mpi::broadcast(comm, info[Controls::status], 0);
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
        LINFO << "Done with PaToH, time : " << MPI_Wtime() - t << "s.";
        t = MPI_Wtime();

        A = CompRow_Mat_double(m_o, n_o, nz_o, val, ir, jc);

        int * test;
        
        nbrows = std::vector<int>(partweights, partweights + icntl[Controls::nbparts]);
        strow =  std::vector<int>(icntl[Controls::nbparts]);

        for(unsigned k = 0; k < icntl[Controls::nbparts]; k++) {
            strow[k] = row_sum;
            row_sum += nbrows[k];
        }

        LINFO << "Finished Partitioning, time: " << MPI_Wtime() - t << "s.";
           
        delete[] ir;
        delete[] jc;
        delete[] val;
        delete[] partvec;
        delete[] partweights;
        delete[] cwghts;
        delete[] pins;
        delete[] xpins;
        delete[] nwghts;
        PaToH_Free();
#else
        info[Controls::status] = -7;
        mpi::broadcast(comm, info[Controls::status], 0);
        throw std::runtime_error("Trying to use PaToH while it is not available");
#endif
        break;
	}
   case 4:
        {
         double t = MPI_Wtime();
         //LINFO << "Reading partition vector " << partvector ;
         int k =  icntl[Controls::nbparts];
         int *partweights = new int[k];

         for(int z =0; z < k; z++){
                partweights[z]=0;
         }

         for(int z =0; z < m_o; z++) {
               partweights[partvec[z]]++;
         }
         row_perm = sort_indexes(partvec, m_o);

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
        nbrows = std::vector<int>(partweights, partweights + icntl[Controls::nbparts]);
        strow =  std::vector<int>(icntl[Controls::nbparts]);

        for(unsigned k = 0; k < icntl[Controls::nbparts]; k++) {
            strow[k] = row_sum;
            row_sum += nbrows[k];
        }
        LINFO << "Done with Manual Partitioning, time : " << MPI_Wtime() - t << "s.";
        break;
        }

    }

    if(write_problem.length() != 0) {
      LINFO << "Writing the problem to the file: " << write_problem;
      int *ir = A.rowptr_ptr();
      int *jc = A.colind_ptr();
      double *val = A.val_ptr();

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

      string parts = write_problem + "_parts";
      f.open(parts.c_str());

      for(unsigned int k = 0; k < (unsigned int)icntl[Controls::nbparts]; k++) {
        f << nbrows[k] << "\n";
      }

      f.close();
    }
}

void abcd::analyseFrame()
{
    LINFO << "Launching frame analysis";
    std::vector<CompCol_Mat_double > loc_parts;
    std::vector<int> ci_sizes;

    double t  = MPI_Wtime();

    column_index.resize(icntl[Controls::nbparts]);
    loc_parts.resize(icntl[Controls::nbparts]);

    LINFO << "Creating partitions";
    
    for (unsigned int k = 0; k < (unsigned int)icntl[Controls::nbparts]; ++k) {
        CompCol_Mat_double part(CSC_middleRows(A, strow[k], nbrows[k]));
        
        int *col_ptr = part.colptr_ptr();
        column_index[k] =  getColumnIndex(col_ptr, part.dim(1));

        // if no augmentation, then create the parts
        if(icntl[Controls::aug_type] == 0)
        {
            parts[k] = CompRow_Mat_double(sub_matrix(part, column_index[k]));
        } else 
        {
            loc_parts[k] = CompCol_Mat_double(part);
        }
    }
    LINFO << "Partitions created in: " << MPI_Wtime() - t << "s.";
    //
#ifdef WIP
    // test augmentation!
    if(icntl[Controls::aug_analysis] == 2){
        double f = 0;
        size_c = 1;
        while(size_c > 0 && f < 0.9){
            dcntl[Controls::aug_filter] = f;
            LDEBUG << "filter value:\t" << fixed << setprecision(5) << f << " gives : ";
            abcd::augmentMatrix(loc_parts);
            f+=0.025;
        }
        exit(0);
    }
#endif //WIP

    if (icntl[Controls::aug_type] != 0) {
        t = MPI_Wtime();
        abcd::augmentMatrix(loc_parts);
        LINFO << "Augmentation time: " << MPI_Wtime() - t << "s.";

        column_index.clear();
        column_index.resize(icntl[Controls::nbparts]);
        ci_sizes.resize(icntl[Controls::nbparts]);
        
        for (unsigned int k = 0;
             k < (unsigned int)icntl[Controls::nbparts]; k++) {

            CompCol_Mat_double &part = loc_parts[k];

            // Build the column index of part
            column_index[k] = getColumnIndex(
                part.colptr_ptr(), part.dim(1)
                );

            std::vector<int> &ci = column_index[k];
            
            ci_sizes[k] = ci.size();

            parts[k] = CompRow_Mat_double(sub_matrix(part, ci));
        }
        if (icntl[Controls::aug_type] != 0)
            LINFO << "Time to regenerate partitions:\t" << MPI_Wtime() - t;
        if (size_c == 0) {
            LWARNING << "WARNING: Size of C is zero, switching to classical cg";
            icntl[Controls::aug_type] = 0;
        }
    }

#ifdef WIP
    // print only the size of C
    if(icntl[Controls::aug_analysis] == 1) exit(0);
#endif // WIP

}

