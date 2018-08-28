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
 * \file utils.cpp
 * \brief Implementation of some usefull utils functions
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include <abcd.h>
#include "blas.h"
#include "mat_utils.h"
#ifndef NO_METIS
    #include <metis.h>
#endif

using namespace std;
using namespace boost::lambda;

/*!
 *  \brief Distribute partitions to masters
 *
 *  Distribution of the partitions to masters:
 *   - #masters == # partitions: 1partition-1master
 *   - 1masters: all partitions for him
 *   - #masters < # partitions:
 *       - if minCommWeight !=0: METIS partitioning of the adjacency graph between partitions
 *      in order to both diminish interconnections and communications. minCommWeight is an
 *      imbalance parameter.
 *       - else: sort partitions from biggest to smallest and distribute to master with current
 *      lowest weight.
 *
 *  \param partitionsSets: groups of partitions for each master
 *  \param weights: vector of weights (#rows) for each matrix
 *  \param nb_parts: number of masters
 *
 */
void abcd::partitionWeights(std::vector<std::vector<int> > &partitionsSets, std::vector<int> weights, int nb_parts)
{
    std::vector<int> sets(nb_parts); // total weight (#rows) of partitions group on each master
    std::map<int, std::vector<int> > pts; // groups of partitions on each master

    for(int i = 0; i < nb_parts; i++){
        sets[i] = 0;
        pts[i];
    }

    // If #masters=#partitions, 1partition per Master
    if (nb_parts == (int)weights.size() ) {
        for(int i = 0; i < nb_parts; i++) {
            pts[i].push_back(i);
            sets[i]+=weights[i];
        }
    // If 1master, it gets all partitions
    } else if (nb_parts == 1) {
        for(int i = 0; i < weights.size(); i++) {
            pts[0].push_back(i);
            sets[0]+=weights[i];
        }
    }
    // If #masters < #partitions
#ifndef NO_METIS
    // minCommWeight: try to diminish comm with METIS or not in config file
    else if (icntl[Controls::minCommWeight] == 0) {
#else
    else {
#endif
	// sort partitions by weights in descending order
	std::vector<int>  sorted = sort_indexes(&weights[0],  weights.size());
        std::sort(weights.begin(), weights.end());
        std::reverse(sorted.begin(), sorted.end());
        std::reverse(weights.begin(), weights.end());
        // distribute the next partition to the current minimum weight master
	for(int i = 0; i < weights.size(); i++){
		int min_index = min_element_index(sets.begin(), sets.end());
		pts[min_index].push_back(sorted[i]);
		sets[min_index] +=  weights[i];
	}
    }
#ifndef NO_METIS
else {
        /* Allocate memory */
        int numparts = weights.size();
        int Wnzmax= numparts*10;
        idx_t Wanz =0;
        // W is an adjacency matrix between partitions
        idx_t *WCp,*WCj;    // Arrays for pointer (p) to columns (j)
        WCp = (idx_t*) calloc(numparts+1,sizeof(idx_t));
        WCj = (idx_t*) calloc(Wnzmax,sizeof(idx_t));
        idx_t *WCx = (idx_t*) calloc(Wnzmax,sizeof(idx_t)); // Array of values
        idx_t *vtxWght = (idx_t*) calloc(numparts,sizeof(idx_t)); // weights from ABCD (#rows in part)

        double timew2=MPI_Wtime();
        idx_t *xb = (idx_t*)calloc(numparts,sizeof(idx_t));
        double *xx = (double*)calloc(numparts,sizeof(double));
        int ind =0;

        /* Build Adjacency matrix between partitions */
        // Find actual columns of the partition
        int ** arrpart = new int*[numparts];
        for(int i=0;i<numparts;i++){
            arrpart[i] = (int*) calloc(A.dim(0), sizeof(int));
            for(int j=0; j<column_index[i].size() ; j++){
                arrpart[i][column_index[i][j]]=1;
            }
        }

        for(int i= 0; i< numparts ; i++){
            vtxWght[i] = weights[i];

            // If array for interaction too small, try to allocate more
            if( (Wanz + numparts) > Wnzmax){
                Wnzmax += numparts*10;
                WCj = (idx_t*)realloc(WCj, sizeof(idx_t)*Wnzmax);
                WCx = (idx_t*)realloc(WCx, sizeof(idx_t)*Wnzmax);
                LINFO << "Array size increased to "<< Wnzmax;
                if(WCj == NULL){cout << "out of Memory" << Wnzmax << " "<< Wanz <<endl;exit(0);}
            }
            // Find interactions
            WCp[i] = Wanz;
            for(int j=0; j< A.dim(0); j++){
                // if part has column, check interactions
                if( arrpart[i][j] ){
                    for(int m=0; m<numparts; m++){
                        if(m == i) continue;
                        // found interaction
                        if( arrpart[m][j] ){
                            // first interaction
                            if(xb[m] != i+1 ){
                                // Add in CSR format
                               xb[m] = i+1; // place of interaction (pointer to col)
                                  WCj[Wanz++] = (idx_t) m;
                                  xx[m] = 1; // # interactions
                            }
                             else{
                                    xx[m]++;
                            }
                        }
                    }
                }
            }
            // Update number of interactions
            for (int pp = WCp[i]; pp < Wanz; pp++){
                WCx[pp] = (idx_t)xx[WCj[pp]];
            }
        }
        WCp[numparts] = Wanz;   // Number of non-zeros for last element in CSR

        // Now you have the SPARSE adjacency matrix

        // Free memory
        for(int i=0;i<numparts;i++)  free(arrpart[i]);
        free(arrpart);

        LINFO << "Preprocessing time for Min. communication partitions distribution " << MPI_Wtime() - timew2;

        /* Build the Graph/adjacency matrix between partitions (interactions/interconnections):
            vertex weights: #rows
            edges weights: #interconnections <=> columns k s.t. there exists lines i,j s.t. aik and ajk non-zeros*/
        // METIS parameters
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_UFACTOR] = 10*icntl[Controls::minCommWeight]; // imbalance ratio (100=10%)
//        options[METIS_OPTION_SEED]    = 1;
//        options[METIS_OPTION_NITER]   = 10;
//        options[METIS_OPTION_NCUTS]   = 1;
//        options[METIS_OPTION_DBGLVL]  = 127;

        idx_t numberofparts = nb_parts;
        idx_t nVertices = numparts;
        idx_t nWeights = 1;
        idx_t *adjncy, *part;
        idx_t objval;   // #cuts
        part = new idx_t[nVertices+1]; // groups of part from METIS
        double t2 = MPI_Wtime();
        int ret;

        // METIS partitioning (Kway or recursive) of the graph
        LINFO << "Launching METIS for Min. communication, aim: " << numberofparts << ", nVertices " << nVertices << ", nEdges " << Wanz ;
        /* METIS: Multilevel algorithm
            - K-way: partition at the coarsest level (but possible empty partitions)
            - Recursive: divide in 2 at coarsest then divides at each level recursively (may be better in some cases)
        */
        ret = METIS_PartGraphKway(&nVertices, &nWeights, WCp, WCj,  vtxWght , NULL, WCx, &numberofparts, NULL,NULL, options, &objval, part);
        //ret = METIS_PartGraphRecursive(&nVertices, &nWeights, WCp, WCj,  NULL, NULL, WCx, &numberofparts, NULL,NULL, options, &objval, part); // same but different heuristic

        LINFO << "Done with METIS, time " << setprecision(6) << MPI_Wtime()-t2 << ", #cuts " << objval << " return " << ret;

        for(int z =0; z < nVertices; z++) {
            pts[part[z]].push_back(z);  // groups of parts in ABCD
            sets[part[z]]+=vtxWght[z];
        }

    }
#endif

    // Sort sets per accumulated weight via tmp_ord to ord_parts
    std::vector<std::pair<int, int>> tmp_ord;
    for (int iii=0; iii<sets.size(); ++iii) {
        std::pair<int, int> p (sets[iii], iii);
        tmp_ord.push_back(p);
    }
    std::sort(tmp_ord.begin(), tmp_ord.end(),
        pair_comparison<int, int,
        position_in_pair::first,
        comparison_direction::descending>);
    std::vector<int> ord_parts;
    for (int iii=0; iii<tmp_ord.size(); ++iii) {
        ord_parts.push_back(tmp_ord[iii].second);
    }

    // Save the sets in ABCD structure
    for(int i = 0; i < nb_parts; i++){
        partitionsSets.push_back(pts[ord_parts[i]]);
    }
}               /* -----  end of function partitionWeights  ----- */

/*!
 *  \brief Parallel dot product p.ap'
 *
 *  Compute the dot product in parallel with interconnected masters thanks to the comm_map
 *
 *  \param p: first vector
 *  \param ap: second vector
 *
 */
double abcd::ddot(VECTOR_double &p, VECTOR_double &ap)
{
    // check corresponding sizes
    int lm = p.size();
    int rm = ap.size();
    if(lm != rm) throw - 800;

    //initialize local values of the product
    VECTOR_double loc_p(lm, 0);
    VECTOR_double loc_ap(rm, 0);
    double loc_r = 0, r = 0;

    int pos = 0;
    for(int i = 0; i < lm; i++) {
        if(comm_map[i] == 1) {
            loc_p(pos) = p(i);
            loc_ap(pos) = ap(i);
            pos++;
        }
    }

    // initialize local value of the product: r=a.ap'
    for(int i = 0; i < lm; i++){
        loc_r += loc_p(i) * loc_ap(i);
    }

    // reduce the value with interconnected processes
    mpi::all_reduce(inter_comm, loc_r, r, std::plus<double>());

    return r;
}               /* -----  end of function abcd::ddot  ----- */

/*!
 *  \brief Compute norm of the current iterate x and the residual b-Ax
 *
 *  Compute norm of the current iterate x and the residual b-Ax
 *
 *  \param x: current iterate
 *  \param b: right hand side
 *  \param nrmR: norm of the current residual
 *  \param nrmX: norm of the current iterate
 *
 */
void abcd::get_nrmres(MV_ColMat_double &x, MV_ColMat_double &b,
                      VECTOR_double &nrmR, VECTOR_double &nrmX)
{
    //initialize local x,r and norms
    int rn = x.dim(1);
    int rm = x.dim(0);

    nrmX = 0;
    nrmR = 0;

    VECTOR_double nrmXV(rn, 0);
    VECTOR_double nrmRV(rn, 0);

    MV_ColMat_double loc_x(rm, rn, 0);
    MV_ColMat_double loc_r(m, rn, 0);
    MV_ColMat_double loc_xfmx(rm, rn, 0);

    // compute local 1-norm of X
    int pos = 0;
    for(int i = 0; i < rm; i++) {
        if(comm_map[i] == 1) {
            for(int j = 0; j < rn; j++) {
                nrmXV(j) += abs(x(i, j));
            }
            pos++;
        }
    }

    pos = 0;
    /* Compute residual loc_r=b-A*x for A on all local partitions */
    for(int p = 0; p < nb_local_parts; p++) {
        for(int j = 0; j < rn; j++) {
            // Compress x
            VECTOR_double compressed_x = VECTOR_double((partitions[p].dim(1)), 0);
            int x_pos = 0;
            for(size_t i = 0; i < local_column_index[p].size(); i++) {
                int ci = local_column_index[p][i];
                compressed_x(x_pos) = x(ci, j);

                x_pos++;
            }
            // Compute loc_r=Ai*x
            VECTOR_double vj =  partitions[p] * compressed_x;
            int c = 0;
            for(int i = pos; i < pos + partitions[p].dim(0); i++)
                loc_r(i, j) = vj[c++];
        }

        pos += partitions[p].dim(0);
    }
    // loc_r=b-Ai*x
    loc_r  = b - loc_r;

    // Compute Inf-norm of the residual
    for(int j = 0; j < rn; j++){
        VECTOR_double loc_r_j = loc_r(j);
        nrmRV(j) = infNorm(loc_r_j);
    }

    // Reduce norms of X and R among all masters
    mpi::all_reduce(inter_comm, nrmRV.ptr(), rn, nrmR.ptr(), mpi::maximum<double>());
    mpi::all_reduce(inter_comm, nrmXV.ptr(), rn, nrmX.ptr(), std::plus<double>());
}               /* -----  end of function abcd::get_nrmres  ----- */

/*!
 *  \brief Binary or
 *
 *  Binary or: return a if non-zero, else b if non-zero, else return zero
 *
 *  \param a: first number
 *  \param b: second number
 *
 */
double or_bin(double &a, double &b){
    if(a!=0) return a;
    else if(b!=0) return b;
    else return 0;
}               /* -----  end of function abcd::or_bin  ----- */

/*!
 *  \brief Sort indexes based on an array and put it in a vector
 *
 *  Sort indexes based on an array and put it in a vector
 *
 *  \param v: array to sort
 *  \param nb_el: size of v
 *
 */
std::vector<int> sort_indexes(const int *v, const int nb_el) {
    // build vector of pairs value-index
    std::vector< std::pair<int,int> > vp;
    vp.reserve(nb_el);
    for(int i = 0; i < nb_el; i++)
        vp.push_back( std::make_pair(v[i], i) );

    // sort indexes based on comparing values in v
    std::sort(vp.begin(), vp.end(),
        pair_comparison<int, int,
        position_in_pair::first,
        comparison_direction::ascending>);
    std::vector<int> idx(vp.size());
    for(int i=0; i<vp.size(); ++i){
      idx[i]=vp[i].second;
    }
    return idx;
}               /* -----  end of function abcd::sort_indexes  ----- */

/*!
 *  \brief Regroups the data from the different sources to a single destination on the root
 *
 *  Centralize the solution vector on the root and compute the forward error if a
 *  starting vector is used for BCG (B=A*Xf with Xf the exact solution).
 *
 *  \param dest: centralized solution vector
 *  \param dest_lda: ???
 *  \param dest_ncols: number of columns in the solution
 *  \param src: local solution vector
 *  \param src_lda: ???
 *  \param src_ncols: number of columns in the local solution
 *  \param globalIndex: correspondance between global and local indices
 *  \param scale: scaling vector on columns
 *
 */
void abcd::centralizeVector(double *dest, int dest_lda, int dest_ncols,
                            double *src,  int src_lda,  int src_ncols,
                            std::vector<int> globalIndex, double *scale)
{
    // Check size of solution in source and destination
    if (src_ncols != dest_ncols) {
        throw std::runtime_error("Source's number of columns must be the same as the destination's");
    }

    MV_ColMat_double source(src, src_lda, src_ncols, MV_Matrix_::ref);

    // Gather local parts of the solution vector on root
    if(IRANK == 0) {
        double t = MPI_Wtime();

        // Centralized solution
        // solution matrix including the output solution being the pointer dest
        MV_ColMat_double vdest(dest, dest_lda, dest_ncols, MV_Matrix_::ref);
        std::map<int, std::vector<double> > xo; // values of the compressed local parts
        std::map<int, std::vector<int> > io; // column indices to uncompress local parts
        std::map<int, int > lo; // position of the local solution part

        // Receive local parts of the solution from other masters
        for(int k = 1; k < inter_comm.size(); k++){
            inter_comm.recv(k, 71, xo[k]);
            inter_comm.recv(k, 72, io[k]);
            inter_comm.recv(k, 73, lo[k]);
        }

        // gather solution parts
        for(int k = 1; k < inter_comm.size(); ++k){
            for(int j = 0; j < dest_ncols; ++j)
                for(size_t i = 0; i < io[k].size() && io[k][i] < n_o; ++i){
                    int ci = io[k][i];
                    vdest(ci, j) = xo[k][i + j * lo[k]] * scale[ci];
                }
        }

        // add part local to the root
        for(int j = 0; j < dest_ncols; ++j)
            for(size_t i = 0; i < globalIndex.size() && glob_to_local_ind[i] < dest_lda; ++i){
                vdest(globalIndex[i], j) = source(i, j) * dcol_[globalIndex[i]];
            }

        // Compute forward error
        ///@TODO Move this away
        if(Xf.dim(0) != 0) {
            MV_ColMat_double xf =  Xf - vdest;
            double nrmxf =  infNorm(xf);
            dinfo[Controls::forward_error] =  nrmxf/nrmXf;
        }
    // Send local part of the solution
    } else {
        std::vector<double> x(src, src + src_lda * src_ncols);

        inter_comm.send(0, 71, x);
        inter_comm.send(0, 72, glob_to_local_ind);
        inter_comm.send(0, 73, src_lda);
    }
}               /* -----  end of function abcd::centralizeVector  ----- */
