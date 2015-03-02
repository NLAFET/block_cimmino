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

#ifndef ABCD_HXX_
#define ABCD_HXX_

#include "mpi.h"

#include <iostream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <map>

#include <mumps.h>
#include "dmumps_c.h"

#include <boost/mpi.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/progress.hpp>
//#include <boost/range/adaptors.hpp>
//#include <boost/range/algorithm.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lexical_cast.hpp>

#include "easylogging++.h"

/*
 * A small hack to make Sparselib++ work with openmpi
 * It's however, better to use it while compiling by putting:
 * -D"COMPLEX=std::complex<double>"
 */
#define COMPLEX std::complex<double>
#include "compcol_double.h"
#include "comprow_double.h"
#include "coord_double.h"
#include "mvm.h"

#include "splib_utils.h"
#include "defaults.h"

/* Some macros*/
#define IRANK inter_comm.rank()
#define IBARRIER inter_comm.barrier()

#define IFMASTER if(inter_comm.rank() == 0)
#define TIC t = MPI_Wtime()
#define TOC MPI_Wtime() - t

#define inter_rank() inter_comm.rank()
#define inter_barrier() inter_comm.barrier()
#define inter_master() if(inter_comm.rank() == 0)

/// Get some goodies from C++11
#define nullptr 0

#ifndef LINFO_
#define LLOG_(v, l) LOG_IF(icntl[Controls::verbose_level] >= v, l)
#define LINFO1 LLOG_(1, INFO)
#define LINFO2 LLOG_(2, INFO)
#define LINFO3 LLOG_(3, INFO)
#define LDEBUG1 LLOG_(1, DEBUG)
#define LDEBUG2 LLOG_(2, DEBUG)
#define LDEBUG3 LLOG_(3, DEBUG)
#endif // LINFO_

#ifdef LINFO
#undef LINFO
#define LINFO LINFO1
#endif //LINFO

#ifdef LDEBUG
#undef LDEBUG
#define LDEBUG LDEBUG1
#endif // LDEBUG

// disable all logging
#ifdef NOLOGGING
#ifdef LINFO
#undef LINFO
#endif // LINFO
#ifdef LDEBUG
#undef LDEBUG
#endif //LDEBUG
#endif //NOLOGGING

using namespace std;
using namespace boost;
using namespace boost::numeric;
using namespace boost::lambda;
//using namespace boost::adaptors;

class abcd
{

public:
    /***************************************************************************
     * Matrix information
     **************************************************************************/
    /*! The number of rows in the matrix */
    int m;

    /*!  The number of columns in the matrix */
    int n; 

    /*! The number of entries in the matrix */
    int nz;

    /*! The symmetry of the matrix
     */
    bool sym;

    /*! The row indices of size #nz */
    int *irn;

    /*! The column indices of size #nz */
    int *jcn;

    /*! Defines wether it's Fortran-Style (1) or C-Style (0, default) */
    int start_index;

    /*! The entries of the matrix of size #nz */
    double *val;

    /*! The number of right-hand sides to solve, default is 1 */
    int nrhs;

    /*! The right-hand side of size #m * #nrhs */
    double *rhs;

    /*!  The solution vector of size #n * #nrhs */
    double *sol;

    /*!  The gateway function that launches all other options
     *
     * Run an operation identified by the value of job_id, it can be
     * either:
     * - -1, initializes the internal matrix used by the
     *      solver. Prior to this call, the user must provide:
     *    * The information about the matrix. #m, #n, #nz,
     *      #sym, #irn, #jcn and #val have to be initialized before the
     *      call.
     *    * After the call, the arrays #irn, #jcn and #val
     *      are no longer used by the solver and can be freely deallocated.
     * - 1, performs the preprocessing. During this call, the solver
     *    scales the matrix, partitions it and, if requested by the user,
     *    performs the augmentation of the matrix. Prior to this call, the
     *    user must provide:
     *  
     *    * The number of partitions to create (see ``nbparts`` )
     *       or ask the solver to guess the appropriate number of
     *       partitions (see ``part_guess``)
     *    * The type of scaling to perform (see ``scaling``)
     *    * The type of augmentation to perform (see ``aug_type``)
     * - 2, creates the augmented systems, analyses them, creates
     *     the mapping between the different mpi-processes and
     *     factorizes the augmented systems.
     * - 3, performs the solution step, the right-hand sides and
     *     their number are required prior to this call.
     *     * The right-hand sides have to be given through the array
     *        ``rhs[]`` and their number in #nrhs.
     *     * The block-size to be used during the block-CG
     *        acceleration. Its value is used only during the regular
     *        block Cimmino solve, and by default its value is 1.
     *     * The solution is centralized (on the master) in the array
     *        ``sol[]``.
     * - 4, combines the call to the phases 1 and 2. 
     * - 5, combines the call to the phases 2 and 3. 
     * - 6, combines the call to the phases 1, 2 and 3.
     *
     */
    int operator() (int job_id);
    
    /**************************************************************************
     * Communication info
    **************************************************************************/
    /*! The integer control array, see Controls::icontrols for the
     *  possible values
     */
    std::vector<int> icntl;

    /*! The real control array, see Controls::dcontrols for the
     *  possible values
     */
    std::vector<double> dcntl;

    /*! The integer info output array, see Controls::info */
    std::vector<int> info;
    /*! The real info output array, see Controls::dinfo */
    std::vector<double> dinfo;

    /**************************************************************************
     * Write problem and log
    **************************************************************************/
    /*! The path where to write the matrix \f$PD_rAD_cP^T\f$ */
    std::string write_problem;

    /*! The path where to write the matrix \f$S_k\f$ where \f$k\f$ is the mpi-process rank */
    std::string write_s;

    /*! The file where to write logging information */
    std::string log_output;
    
    /**************************************************************************
     * Partitioning informations
    **************************************************************************/
    /*! The starting row index of each partition */
    std::vector<int> strow;

    /*! The number of rows per partition */
    std::vector<int> nbrows;

    /**************************************************************************
     * Communication info
    **************************************************************************/
    /*! The global communicator */
    mpi::communicator comm; 
    /*! The number of parallel CG instances */
    int parallel_cg;

    int initializeMatrix();
    int preprocessMatrix();
    int factorizeAugmentedSystems();
    int solveSystem();
    
    abcd();
    ~abcd();

    // 1-based vectors representing local S matrix in coordinate format
    std::vector<int> S_rows;
    std::vector<int> S_cols;
    std::vector<double> S_vals;
    // shape of S is the size of C augmentation part
    int size_c;

    // the residual info
    std::vector<double> rhoVector;
    std::vector<double> scaledResidualVector;

private:
    int gqr(MV_ColMat_double &P, MV_ColMat_double &AP, MV_ColMat_double &R, int s, bool use_a);
    int gqr(MV_ColMat_double &p, MV_ColMat_double &ap, MV_ColMat_double &r, CompCol_Mat_double g, int s, bool use_a);
    void gmgs2(MV_ColMat_double &p, MV_ColMat_double &ap, MV_ColMat_double &r, int s, bool use_a);
    void gmgs2(MV_ColMat_double &p, MV_ColMat_double &ap, MV_ColMat_double &r, CompCol_Mat_double g, int s, bool use_a);

    // norms
    double nrmA;
    std::vector<double> nrmB;
    double nrmXf;
    double nrmMtx;

    // to check the order of calls
    int last_called_job;

    // preprocess stuffs
    std::vector<double> drow_;
    std::vector<double> dcol_;

    void scaling();

    // Scales the matrix
    void scaleMatrix(int norm);
    void diagScaleMatrix(std::vector<double> & , std::vector<double> & );
    void diagScaleRhs(VECTOR_double &);
    void diagScaleRhs(MV_ColMat_double &);

    // structure functions
    /// Partitions the matrix into abcd::nbrows
    void partitionMatrix();
    /**
     * Analyses the structure of each partition
     * Compresses the  and analyses the interconnections between them
     */
    void analyseFrame();

    /*! The type of process
     * - 0 the process will behave as a CG master
     * - 1 the process will be a  MUMPS Slave
     */
    int instance_type;
    
    // Build the augmented version of the matrix (ABCD)
    void augmentMatrix(std::vector<CompCol_Mat_double > &loc_parts);
    void cijAugmentMatrix(std::vector<CompCol_Mat_double > &loc_parts);
    void aijAugmentMatrix(std::vector<CompCol_Mat_double > &loc_parts);

    // Communication stuffs
    void createInterCommunicators();
    void distributePartitions();
    void createInterconnections();

    void distributeData();

    void solveABCD(MV_ColMat_double &b);
    MV_ColMat_double solveS ( MV_ColMat_double &f );

    void buildS(std::vector<int> &rows,
                std::vector<int> &cols,
                std::vector<double> &vals);

    void buildS(std::vector<int> &rows,
                std::vector<int> &cols,
                std::vector<double> &vals,
                std::vector<int> &columns_to_build);

    Coord_Mat_double buildS();
    Coord_Mat_double buildS(std::vector<int>);

    // Cimmino
    void initializeDirectSolver();
    void distributeRhs();
    void distributeNewRhs();
    void bcg(MV_ColMat_double &b);

    MUMPS buildM();
    VECTOR_double solveM ( MUMPS &mu, VECTOR_double &z );
    MV_ColMat_double prodSv(MV_ColMat_double &);
    VECTOR_double pcgS ( VECTOR_double &b );
    std::vector<int> selected_S_columns;
    std::vector<int> skipped_S_columns;
    double compute_rho(MV_ColMat_double &X, MV_ColMat_double &U);

    // MUMPS
    int m_n;
    int m_nz;
    int n_aug, nz_aug;
    std::vector<int> irn_aug, jcn_aug;
    std::vector<double> val_aug;

    MUMPS mumps;
    void initializeMumps(MUMPS &, bool local);
    void initializeMumps(MUMPS &);
    void createAugmentedSystems(int &n_aug,
                                int &nz_aug,
                                std::vector<int> &irn_aug,
                                std::vector<int> &jcn_aug,
                                std::vector<double> &val_aug);
    void analyseAugmentedSystems(MUMPS &);
    void allocateMumpsSlaves(MUMPS &);
    void factorizeAugmentedSystems(MUMPS &);

    MV_ColMat_double sumProject(double alpha,
                                MV_ColMat_double &Rhs,
                                double beta,
                                MV_ColMat_double &X);
    MV_ColMat_double spSimpleProject(std::vector<int> mycols);
    void spSimpleProject(std::vector<int> mycols, std::vector<int> &vrows,
                         std::vector<int> &vcols, std::vector<double> &vvals);

    int my_master;
    std::vector<int> my_slaves;


    void waitForSolve();
    std::vector<int> comm_map;

    // SOme utilities
    void partitionWeights(std::vector<std::vector<int> > &parts,
                          std::vector<int> weights, int nb_parts);
    void partitioning(std::vector<std::vector<int> > &, std::vector<int>, int);
    double ddot(VECTOR_double &p, VECTOR_double &ap);
    void get_nrmres(MV_ColMat_double &x,
                    MV_ColMat_double &b,
                    VECTOR_double &nrmR,
                    VECTOR_double &nrmX);

    MUMPS mumps_S;
    Coord_Mat_double S;
    inline int S_nnz() { return S_vals.size(); }

    std::map<int, CompRow_Mat_double> parts;
    std::vector<CompRow_Mat_double> partitions;
    std::vector<std::vector<int> > partitionsSets;

    MV_ColMat_double Xf;
    MV_ColMat_double B;
    MV_ColMat_double RRand;

    CompRow_Mat_double A;
    std::vector<int> row_perm;

    bool runSolveS;

    int n_l, m_l, nz_l;
    int n_o, m_o, nz_o;

    MV_ColMat_double Xk;

    int nbparts; /// The number of partitions
    int nb_local_parts;

    /// A reverse index of columns, contains the original index of each column for each partition
    std::vector<std::vector<int> > column_index;
    std::vector<std::map<int,int> > column_index_cache;
    /// A merge of col_index vectors, determines non-null columns in all local partitions
    std::vector<std::vector<int> > local_column_index;
    int **fast_local_column_index;
    std::map<int,int> glob_to_local;
    // Both of these are meant to replace the glob_to_local map
    std::vector<int> glob_to_local_ind;
    std::map<int,int> glob_to_local_c;
    std::vector<int>::iterator st_c_part_it;
    int st_c_part;

    std::vector<std::map<int,int> > glob_to_part;
    std::vector<std::map<int,int> > part_to_glob;
    std::vector<int> stC;
    bool use_xk;
    bool use_xf;
    MV_ColMat_double solution; ///< The solution vector

    /**
     * Contains the mutual interconnections between partitions
     * The key is the cg-master rank (in inter_comm) and the value is the column indices
     */
    std::map<int,std::vector<int> > col_interconnections;
    std::map<int,std::vector<int> > col_inter;
    /// Contains the partitions that are handled by this instance
    std::vector<int> parts_id;

    // logging stuffs
    bool verbose;
    // easyloggingpp::Configurations log_config;

    /// The communicator shared by CG masters
    mpi::communicator inter_comm;
    /// The communicator of local slaves
    mpi::communicator intra_comm; 

    void centralizeVector(double *dest, int dest_lda, int dest_ncols,
                          double *src,  int src_lda,  int src_ncols,
                          std::vector<int> globalIndex, double *scale);

};


void configure_logger(std::string log_file);
void logger_set_filename(std::string log_file);
    
double or_bin(double &a, double &b);
std::vector<int> sort_indexes(const int *v, const int nb_el);

#endif // ABCD_HXX
