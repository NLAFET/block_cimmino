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
 * \file structure.cpp
 * \brief Implementation of the functions computing the partitions
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include <abcd.h>
#include "mat_utils.h"
#include <iostream>
#include <fstream>
#include <vect_utils.h>
#include <boost/lambda/lambda.hpp>

#ifdef PATOH
#include <patoh.h>
#endif //PATOH

#include <metis.h>

using namespace boost::lambda;

/*!
 *  \brief Partition the Matrix
 *
 *  Guess the number of partitions if needed and partition the matrix using the chosen
 *  method:
 *   - 1: uniform partitioning with given nbrows
 *   - 2: uniform partitioning
 *   - 3: PaToH partitioning
 *   - 4: manual partitioning input
 *  Then compute overlapping and write the scaled and permuted matrix.
 *
 */

CompRow_Mat_double ColumnSparsifying(CompRow_Mat_double &inMtx){

	// Parameter for dropping threshold
    int metis_ColSparsing = 1;

    if(metis_ColSparsing <= 0)
        return inMtx;

    int m,i,j,k, *Ap, *Aj ;
    double *Ax;
    long long newVal, fre;

    m = inMtx.dim(0);
    Ap = inMtx.rowptr_ptr();
    Aj = inMtx.colind_ptr();
    Ax = inMtx.val_ptr();

    int scaleMax= 1000; // number of bins
    int *freqs = new int [scaleMax+1]();

    int nnz_A = Ap[m];
    int droppedRow=0;
    int droppedCol=0;

    int *Cp = new int [m + 1];
    int *Cj = new int[nnz_A];
    double *Cx = new double[nnz_A];
    CompRow_Mat_double Amtx;

    for(i=0;i<m+1;i++){
        Cp[i] = 0;
    }

    //Column Sparsing
    if(metis_ColSparsing > 0){
        double sq = sqrt(m) * metis_ColSparsing;
        Amtx = csr_transpose(inMtx); // make it column oriented
        Ap = Amtx.rowptr_ptr();
        Aj = Amtx.colind_ptr();
        Ax = Amtx.val_ptr();

        for(i=0;i<m;i++){
            if( (Ap[i+1] - Ap[i] ) > sq)
            {
                droppedCol++;
                int todroppedNNZ = (Ap[i+1] - Ap[i]) - sq;

                double max = -DBL_MAX;
				double min = DBL_MAX;
                
				// finding min and max in i^th column
				for(j=Ap[i];j<Ap[i+1];j++){
                    if(abs(Ax[j]) < min && Ax[j] != 0) {
                        min = abs(Ax[j]);
                    }
                    if (abs(Ax[j]) > max){
                        max = abs(Ax[j]);
                    }
                }

		// linear distribution for values
                double scaleRatio = abs((scaleMax) / (max-min)); 

		// initializing
                for(int ff=0;ff<scaleMax;ff++) freqs[ff]=0;

		// find frequencies of each bin
                Cp[i+1] = Cp[i];
                for(j=Ap[i];j<Ap[i+1];j++){
                    newVal = floor(scaleRatio * abs((abs(Ax[j])-min)));
                    freqs[newVal]++;
                }

		// deciding whichh bins will be dropped 
                int sum=0;
                for(fre=0; fre < scaleMax;fre++){
                    sum += freqs[fre];
                    if(sum > todroppedNNZ )
                        break;
                }
		
		// drop values in bins < "fre" until "todroppedNNZ=0"
                for(j=Ap[i];j<Ap[i+1];j++){
                    newVal = floor(scaleRatio * abs((abs(Ax[j])-min)));

                    if((newVal > fre) || todroppedNNZ == 0 ){
                        k = Cp[i+1]++;
                        Cj[k] = Aj[j];
                        Cx[k] = Ax[j];
                    }
                    else{
                        todroppedNNZ--;
                    }
                }
             }
             else{
		// ToDo: this can be improved
		// for already sparse columns, copy original column without sparsifying
                Cp[i+1] = Cp[i];
                for(j=Ap[i];j<Ap[i+1];j++){
                    k = Cp[i+1]++;
                    Cj[k] = Aj[j];
                    Cx[k] = Ax[j];
                }
            }
        }
    }

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    if(myid==0) cout << "Column Dropped Cp[m] " << Cp[m] << " DroppedCol "<< droppedCol << " " <<  "metis_ColSparsing factor:"<< metis_ColSparsing  << endl;

    CompRow_Mat_double E (m , m, Cp[m] , Cx, Cp, Cj);

    delete Cp, Cj,Cx;
    delete freqs;

    if(metis_ColSparsing > 0){
           return  csr_transpose(E); //CSC to CSR
     }
    return E;
}




void abcd::partitionMatrix()
{
    unsigned handled_rows = 0;
    unsigned ceil_per_part, floor_per_part;
    unsigned row_sum = 0;
    // If the user wants #parts to be guessed then:
    // m-#parts is: [1]-1; [2-8]-2; [9-1 000]-4; [1001-50 000]-8;
    //              then [<100 000]-10 000 per part
    //              else 20 000 per part
    int guessPartitionsNumber = icntl[Controls::part_guess];
    if(guessPartitionsNumber == 1 && icntl[Controls::part_type] > 1 && icntl[Controls::part_type] != 4){
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

    /* CHECKS */
    // check #parts
    if (icntl[Controls::nbparts] <= 0){
        info[Controls::status] = -3;
        mpi::broadcast(comm, info[Controls::status], 0);
        throw std::runtime_error("FATAL ERROR: Number of partitions is negative or zero");
    }
    // check #masters <= #parts
    if (icntl[Controls::nbparts] < parallel_cg) {
        // the user is not supposed to know about the parallel_cg
        // LERROR << "ERROR: Number of partitions is smaller than the number of parallel_cg";
        LERROR << "Oops! This should not happen";
        LWARNING << "WARNING: Increasing the number of partitions from " << icntl[Controls::nbparts]
                 << " up to " << parallel_cg;
        icntl[Controls::nbparts] = parallel_cg;
    }
    // check #parts <= m
    if(icntl[Controls::nbparts] > m) {
        LERROR << "ERROR: Number of partitions is larger than the number of rows";
        LWARNING << "WARNING: Decreasing the number of partitions from " << icntl[Controls::nbparts]
                 << " down to " << m;
        icntl[Controls::nbparts] = m;
    }
    // if 1 part then no PaToH
    if (icntl[Controls::nbparts] == 1 && icntl[Controls::part_type] == 3) {
        LWARNING << "WARNING: PaToH is useless with a single partiton request, switching to automatic partitioning";
        icntl[Controls::part_type] = 2;
    }

    /* Perform specified partitioning */
    switch(icntl[Controls::part_type]){
        /*-----------------------------------------------------------------------------
         *  Uniform partitioning with a given nbrows
         *-----------------------------------------------------------------------------*/
    case 1: {
        int row_sum=0;
        for(unsigned k = 0; k < (unsigned)icntl[Controls::nbparts]; k++) {
      	    row_indices.push_back(vector<int>());
            if (nbrows[k] == 0) {
                info[Controls::status] = -6;
                mpi::broadcast(comm, info[Controls::status], 0);
                throw std::runtime_error("FATAL ERROR: You requested an empty partition.");
            }
	    for(int ii=row_sum; ii<row_sum + nbrows[k];ii++ ){
		row_indices[k].push_back(ii);
	    }
            row_sum += nbrows[k];
        }
        break;
        /*-----------------------------------------------------------------------------
         *  Uniform partitioning with only icntl[Controls::nbparts] as input (generates nbrows)
         *-----------------------------------------------------------------------------*/
    } case 2: {
        unsigned floor_per_part;
        floor_per_part = floor(float(m_o)/float(icntl[Controls::nbparts]));
        int remain = m_o - ( floor_per_part * icntl[Controls::nbparts]);
        nbrows = std::vector<int>(icntl[Controls::nbparts]);
        for(unsigned k = 0; k < (unsigned) icntl[Controls::nbparts]; k++) {
           int cnt = floor_per_part;
           if(remain >0){
               remain--;
               cnt++;
           }
           nbrows[k] = cnt;
        }
        int row_sum=0;
        for(unsigned k = 0; k < (unsigned)icntl[Controls::nbparts]; k++) {
       	    row_indices.push_back(vector<int>());
	    for(int ii=row_sum; ii<row_sum + nbrows[k];ii++ ){
		row_indices[k].push_back(ii);
	    }
            row_sum += nbrows[k];
        }
        break;
        /*-----------------------------------------------------------------------------
         *  PaToH partitioning
         *-----------------------------------------------------------------------------*/
    } case 3: {
#ifdef PATOH
	/* Initialize PaToH parameters */
        PaToH_Parameters args;
        int _c, _n, _nconst, _ne, *cwghts, *nwghts, *xpins, *pins, *partvec,
            cut, *partweights, ret;
	double _imba;

        CompCol_Mat_double t_A = Coord_Mat_double(A);

        double t = MPI_Wtime();
        LINFO << "Launching PaToH";

        PaToH_Initialize_Parameters(&args, PATOH_CUTPART, PATOH_SUGPARAM_DEFAULT);
        args._k = icntl[Controls::nbparts];

	/* Allocate PaToH structure */
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

	/* Launch PaToH partitioning */
        args.final_imbal    = _imba;
        args.init_imbal     = _imba * 2.0;
        args.seed           = 1;

        partvec     = new int[_c];
        partweights = new int[args._k * _nconst];

        PaToH_Part(&args, _c, _n, _nconst, 0, cwghts, nwghts,
                   xpins, pins, NULL, partvec, partweights, &cut);

        // Check no empty partitions
	for (int i = 0; i < icntl[Controls::nbparts]; i++) {
	    // row_indices: stores indices of rows for each part
       	    row_indices.push_back(vector<int>());

            if (partweights[i] == 0) {
                info[Controls::status] = -6;
                mpi::broadcast(comm, info[Controls::status], 0);
                throw std::runtime_error("FATAL ERROR: PaToH produced an empty partition, Try to reduce the imbalancing factor");
            }
        }

	for(int zz=0;zz<_c;zz++){
		row_indices[partvec[zz]].push_back(zz);
	}
        LINFO << "Done with PaToH, time : " << MPI_Wtime() - t << " s.";
        t = MPI_Wtime();

        LINFO << "Finished Partitioning, time: " << MPI_Wtime() - t << " s.";

        /* Free memory */
        //delete[] ir;
        //delete[] jc;
        //delete[] val;
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
        /*-----------------------------------------------------------------------------
         *  Manual partitionning input
         *-----------------------------------------------------------------------------*/
   } case 4: {
        double t = MPI_Wtime();
        int k =  icntl[Controls::nbparts];
        int *partweights = new int[k];

        // use partvec to initialize partitions
        for(int z =0; z < k; z++){
	    row_indices.push_back(vector<int>());
            partweights[z]=0;
        }

        for(int z =0; z < m_o; z++) {
            if(partvec[z] > k-1) {
                stringstream err;
                LERROR << "Error : Appears to be too many partitions in manual partitioning file.";
                throw std::runtime_error(err.str());
            }
            partweights[partvec[z]]++;
	    row_indices[partvec[z]].push_back(z);
        }

        // check no empty partitions
        for(int i=0; i<k; i++) {
            if(!partweights[i]) {
                stringstream err;
                LERROR << "Error : A partition from manual partitioning was found empty," <<
                    "check your partitioning file.";
                throw std::runtime_error(err.str());
            }
        }

        LINFO << "Done with Manual Partitioning, time : " << MPI_Wtime() - t << "s.";
        break;
    } 
    case 5: {
        double t = MPI_Wtime();

        CompRow_Mat_double Ddrop = ColumnSparsifying(A);
        LINFO << "Sparsifying Columns Time " << MPI_Wtime() - t;

        CompRow_Mat_double AAT;
        CompRow_Mat_double AT = csr_transpose(Ddrop);

        double tt2 = MPI_Wtime();
        AAT = spmm_overlap(Ddrop,AT);
        double tt3 = MPI_Wtime();

        int nnz_AAT = AAT.NumNonzeros();
        LINFO<< "After SPMM nnz on AAT " <<  nnz_AAT << " SPMM Time " << tt3-tt2 ;

        int i,k,j;
        int *jiro = AAT.rowptr_ptr();
        int *jjco = AAT.colind_ptr();
        double *jvalo = AAT.val_ptr();

        idx_t *val_;
        idx_t *adjncy, *part;
        idx_t objval;
        idx_t *ewgt, *vwgt;
        idx_t *xadj  = new idx_t [n_o + 1];

        val_ = new idx_t[nnz_AAT];

	// for dropping some tiny values in graph after AAT
	double metis_epsilon = 0; // 1e-16;

	// scale factor for casting integer
	int metis_scaleMax = 1000000;

	// for enabling edge and vertex weights
	int metis_vwght = 0; // 0: balance on #rows
	int metis_ewght = 1;

	double scaleRatio = abs(metis_scaleMax) ;

	
	// for casting floating number to integer
        for(int i = 0; i < nnz_AAT; i++){
            if(jvalo[i]==0) val_[i] = 0;
            else{
		// assuming jvalo is taken its absolute value during SPMM
                double newVal = scaleRatio *  jvalo[i] ;
                if(newVal < metis_epsilon){
                    val_[i] = 0;
                }
                else
                    val_[i]  = ceil(newVal);
                }
        }
        
	int *Ap = A.rowptr_ptr();
        if(metis_vwght != 0)
        {
            vwgt = new idx_t[n_o+1];
            for(int i = 0 ; i < n_o; i++){
                vwgt[i] = (int)(sqrt(Ap[i+1] - Ap[i])/2)+1;
                if (vwgt[i] == 0 )  cout << " ERROR vertex weight= 0 " << endl;
            }
        }
        else
        {
            vwgt=NULL;
        }

        for(int i = 0; i <= n_o; i++){
            xadj[i] = 0;
        }

        adjncy = new idx_t[nnz_AAT];
        double *tmp_val = new double[nnz_AAT];

        int nz_new = 0;

	// eliminate zero and diagonal values then fill tmp_val
        for (k = 0; k < n_o; k++) {
            xadj[k+1] = xadj[k];
            for (j =jiro[k]; j < jiro[k+1]; j++) {
                //if(val_[j] != 0 && jjco[j] != k)
                if(val_[j] != 0)
                {
                    i = xadj[k+1]++;
                    adjncy[i] = jjco[j];
                    if( metis_ewght != 0)
                    {
                        tmp_val[i] = val_[j];
                    }
                }
            }
        }
        nz_new = xadj[k];

        if(metis_ewght != 0)
        {
            ewgt = new idx_t[nz_new];
            for(int i = 0; i < nz_new; i++){
                ewgt[i] = tmp_val[i];
            }
        }
	else {
            ewgt = NULL;
        }

        delete[] tmp_val;
        delete[] val_;

        printf("Graph scaleFactor %d graphEpsilon %e METIS Preprocessing Time %lf\n",metis_scaleMax, metis_epsilon,MPI_Wtime()-tt3);

        part = new idx_t[n_o+1];
        idx_t nParts = icntl[Controls::nbparts];

        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
	
		int metis_seed = -1;
		int metis_ufactor = 100; // 10% imbalance
		int metis_recursive = 0; // for calling METIS_PartGraphRecursive (optional)

        options[METIS_OPTION_SEED]    = metis_seed;
        options[METIS_OPTION_UFACTOR] = metis_ufactor;
        //options[METIS_OPTION_DBGLVL]  = 17;

        idx_t nVertices = n_o;
        idx_t nWeights = 1;

        LINFO << "Launching METIS " << " Vwgt: " << metis_vwght <<  " Ewgt: " << metis_ewght  <<" UFAC: " << metis_ufactor << " Recursive: " << metis_recursive << " part: " << nParts;

        double t2 = MPI_Wtime();
        int ret;
 	   
        // for calling METIS_PartGraphRecursive (optional)
        if(metis_recursive == 0){
            ret = METIS_PartGraphKway(&nVertices,&nWeights, xadj, adjncy,  vwgt, NULL, ewgt, &nParts, NULL,NULL, options, &objval, part);
        }else{
            ret = METIS_PartGraphRecursive(&nVertices,&nWeights, xadj, adjncy,  vwgt, NULL, ewgt, &nParts, NULL,NULL, options, &objval, part);
        }

        LINFO << "Done with METIS, time " << setprecision(6) << MPI_Wtime()-t2 << " cut " << objval << " return " << ret;
        int *partweights = new int[nParts];
        int *partvec = new int[n_o];
        for(int z =0; z < nParts; z++)
        {
            partweights[z]=0;
        }
        for(int z =0; z < n_o; z++)
        {
            partweights[part[z]]++;
        }


		// Check no empty partitions
        for (int i = 0; i < icntl[Controls::nbparts]; i++) {
           // row_indices: stores indices of rows for each part
            row_indices.push_back(vector<int>());
            if (partweights[i] == 0) {
                info[Controls::status] = -6;
                mpi::broadcast(comm, info[Controls::status], 0);
                throw std::runtime_error("FATAL ERROR: METIS produced an empty partition, Try to reduce the imbalancing factor");
            }
        }
 
        for(int zz=0;zz<n_o;zz++){
                row_indices[part[zz]].push_back(zz);
        }


         delete[] part;
         delete[] partvec;
         delete[] partweights;
         if(ewgt!=NULL) delete[] ewgt;
         if(vwgt!=NULL) delete[] vwgt;
         delete[] xadj;
         delete[] adjncy;

         LINFO << "Done with Numerically Aware Partitioning, time : " << setprecision(6) << MPI_Wtime() - t << " s.";
	break;
    }
    default: {
        // Wrong job id
        info[Controls::status] = -1;
        throw std::runtime_error("Wrong partitioning type.");
    }
  }

    /* Overlap num_overlap rows at the start and at the end of each partition */
    if  (icntl[Controls::num_overlap] > 0){
        if (icntl[Controls::overlap_strategy] == 1){
            partvec = new int[n_o];
            for(unsigned k = 0; k < icntl[Controls::nbparts]; k++) {
                for(unsigned r = 0; r < row_indices[k].size() ; r++) {
                    partvec[row_indices[k][r]] = k; //row_indices[k][r];
                    //~ cout << r << " " << partvec[row_indices[k][r]] << endl;
                }
            }

            // AT = transpose of A
            CompRow_Mat_double AT = csr_transpose(A);

            // Efficient Sequential Matrix Multiplication AAT = A * A'
            CompRow_Mat_double AAT = spmm_overlap(A,AT);

            int *jiro = AAT.rowptr_ptr();
            int *jjco = AAT.colind_ptr();
            double *jvalo = AAT.val_ptr();
            for(int k = 0; k < icntl[Controls::nbparts]; k++) {
                std::vector<std::pair<int, double>> selected_rep;
                for(int r = 0; r <  n_o  ; r++) {
                    if( partvec[r] != k) {
                        double sum =0;
                        for(int j =jiro[ r ]; j< jiro[ r+1];j++){
                            if( partvec[ jjco[j] ] == k) {
                                sum += jvalo[j];
                            }
                        }
                        selected_rep.push_back(std::make_pair(r,sum));
                    }
                }
                std::sort(selected_rep.begin(), selected_rep.end(),   pair_comparison<int, double,  position_in_pair::second,  comparison_direction::descending>);

                for(int jj = 0; jj<icntl[Controls::num_overlap]; jj++){
                    row_indices[k].push_back ( selected_rep[jj].first );
                    LINFO3 << "Row " << selected_rep[jj].first << " added ("<< selected_rep[jj].second <<") as "<< jj+1 <<". replicated row into " << k << " .part";
                }
            }
	} else if (icntl[Controls::overlap_strategy] == 0){
            LINFO2 << "Duplicating " << icntl[Controls::num_overlap] <<
                " overlapping rows between partitions.";
            for(unsigned k = 0; k < icntl[Controls::nbparts]; k++) {
                for(int zz = 0; zz < icntl[Controls::num_overlap]; zz++) {
                    if(k < icntl[Controls::nbparts]-1 ){
                        if(icntl[Controls::num_overlap] >=  row_indices[k+1].size()){
                            throw std::runtime_error(" More #rows replication than the successive block");
                        } else{
                            row_indices[k].push_back( row_indices[k+1][zz]  );
                            LINFO3 << "Row " << row_indices[k+1][zz] << " added as "<< zz+1 <<". replicated row into " << k << " .part";
                        }
                    } else if( k > 0){
                        if(icntl[Controls::num_overlap] >=  row_indices[k-1].size()){
                            throw std::runtime_error(" More #rows replication than the previous block");
                        } else{
                            row_indices[k].push_back( row_indices[k-1][zz]  );
                            LINFO3 << "Row " << row_indices[k-1][zz] << " added as "<< zz+1 <<". replicated row into " << k << " .part";
                        }
                    }
                }
            }
        }
    }

    /* Write the (permuted) matrix A and the its uniform partitioning (#rows1, ...) to a file */
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
            f << row_indices[k].size() << "\n";
        }

        f.close();
    }
}    /* ----- end of method abcd::partitionMatrix ----- */

/*!
 *  \brief Analyses/Compress partitions
 *
 *  Analyses the structure of each partition, in particular finds the non-empty columns,
 *  then compresses it and finds the interconnected columns with other partitions.
 *  In the case of augmentation, augment the partitions and also find the interconnections.
 *
 */
void abcd::analyseFrame()
{
    LINFO << "Launching frame analysis";

    double t  = MPI_Wtime();

    column_index.resize(icntl[Controls::nbparts]);
    std::vector<CompCol_Mat_double > loc_parts;
    loc_parts.resize(icntl[Controls::nbparts]);

    LINFO << "Creating partitions";

    // Create the compressed partitions
    for (unsigned int k = 0; k < (unsigned int)icntl[Controls::nbparts]; ++k) {
        CompCol_Mat_double part(CSC_extractByIndices(A, row_indices[k] ));

        // Find non-empty columns
        int *col_ptr = part.colptr_ptr();
        column_index[k] =  getColumnIndex(col_ptr, part.dim(1));

        // Create the part (if augmentation, just a local tmp part)
        if(icntl[Controls::aug_type] == 0) {
            parts[k] = CompRow_Mat_double(sub_matrix(part, column_index[k]));
        } else {
            loc_parts[k] = CompCol_Mat_double(part);
        }
    }
    LINFO << "Partitions created in: " << MPI_Wtime() - t << "s.";

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
        std::vector<int> ci_sizes;

        t = MPI_Wtime();
        // Augment the compressed partitions
        abcd::augmentMatrix(loc_parts);
        LINFO << "Augmentation time: " << MPI_Wtime() - t << "s.";

        // Recompute column_index
        column_index.clear();
        column_index.resize(icntl[Controls::nbparts]);
        ci_sizes.resize(icntl[Controls::nbparts]);

        for (unsigned int k = 0;
             k < (unsigned int)icntl[Controls::nbparts]; k++) {
            // Build the column index of the local partition
            CompCol_Mat_double &part = loc_parts[k];

            // Build the column index of part
            column_index[k] = getColumnIndex(
                part.colptr_ptr(), part.dim(1)
                );

            // Save the compressed augmented partitions
            std::vector<int> &ci = column_index[k];
            ci_sizes[k] = ci.size();
            parts[k] = CompRow_Mat_double(sub_matrix(part, ci));
        }

        if (icntl[Controls::aug_type] != 0)
            LINFO << "Time to regenerate partitions:\t" << MPI_Wtime() - t;

        // If no actual augmentation, switch to classic BCG
        if (size_c == 0) {
            LWARNING << "WARNING: Size of C is zero, switching to classical cg";
            icntl[Controls::aug_type] = 0;
        }
    }

#ifdef WIP
    // print only the size of C
    if(icntl[Controls::aug_analysis] == 1) exit(0);
#endif // WIP
}               /* -----  end of function abcd::analyseFrame  ----- */
