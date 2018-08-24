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

// the configured options and settings for ABCD

/*!
 * \file abcd.h
 * \brief Header for the class ABCD
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#define ABCD_VERSION_MAJOR 1
#define ABCD_VERSION_MINOR 0

#ifndef ABCD_HXX_
#define ABCD_HXX_

#include "mpi.h"
#if defined(_OPENMP)
   #include <omp.h>
#endif

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
#include <boost/serialization/utility.hpp>

#include "easylogging++.h"

#include "comparison_op.h"

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
/*******************************************************************************
 * Attributes
 ******************************************************************************/
    /***************************************************************************
     * Matrix information
     **************************************************************************/
    /*! The number of rows in the local matrix part */
    int m;

    /*!  The number of columns in the local matrix part */
    int n;

    /*! The number of entries in the local matrix part */
    int nz;

    /*! The saved number of rows, columns, non-zero values of the original complete matrix */
    int n_o, m_o, nz_o;

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

    /*! The starting point for BCG */
    MV_ColMat_double Xk;

    /* Specify if a starting point for BCG is used */
    bool use_xk;

    /* the residual vector */
    std::vector<double> rhoVector;

    /* the scaled residual vector */
    std::vector<double> scaledResidualVector;

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
    /*! The real scaling Number of iterations, see Controls::scaling */
    std::vector<int> man_scaling;

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
    /*! row indices of partitions for each partition */
    vector< vector<int> > row_indices;

    /*! The number of rows per partition */
    std::vector<int> nbrows;

     // for reading partvector
    int *partvec;

    /**************************************************************************
     * Communication info
    **************************************************************************/
    /*! The global communicator */
    mpi::communicator comm;

    /*! The number of parallel CG instances (#masters) */
    int parallel_cg;

    /**************************************************************************
     * Augmentation
    **************************************************************************/
    /* 1-based vectors representing local S matrix in coordinate format */
    std::vector<int> S_rows;
    std::vector<int> S_cols;
    std::vector<double> S_vals;

    /* shape of S is the size of C augmentation part */
    int size_c;

/*******************************************************************************
 * Methods
 ******************************************************************************/
    abcd();
    ~abcd();
    int initializeMatrix();
    int preprocessMatrix();
    int factorizeAugmentedSystems();
    int solveSystem();
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

private:
/*******************************************************************************
 * Attributes
 ******************************************************************************/
    /**************************************************************************
     * General
    **************************************************************************/
    /* to check order of phases */
    int last_called_job;

    /**************************************************************************
     * Preprocessing
    **************************************************************************/
    /* Matrix (possibly augmented) */
    CompRow_Mat_double A;
//    double nrmA;

    /* dimensions of the (augmented or not) matrix */
    int n_l, m_l, nz_l;

    /* row/column scaling factors to scale the system */
    std::vector<double> drow_;
    std::vector<double> dcol_;

    /* norm of the preprocessed matrix (possibly augmented) */
    double nrmMtx;

    /* Specify if an artificial RHS must be used */
    bool use_xf;

    /* Artificial RHS, which columns should be of values (i+1)/m before scaling */
    /* Then solution vectors for BCG */
    MV_ColMat_double Xf;

    /* norm of the matrix Xf used as basis for the artificial right hand side (B=A*Xf) */
    double nrmXf;

    /* Right Hand Side (local) */
    MV_ColMat_double B;
//    MV_ColMat_double RRand;

    /* Infinite norm of the Right Hand Side */
    std::vector<double> nrmB;

    /* A merge of col_index vectors, determines non-null columns in all local partitions */
    std::vector<std::vector<int> > local_column_index;
//    int **fast_local_column_index;

    /* The solution vector */
    MV_ColMat_double solution;

    /**************************************************************************
     * Partitioning
    **************************************************************************/
    /* The number of partitions */
    int nbparts;

    /* Map containing for each partition the index of its non-empty columns */
    std::vector<std::vector<int> > column_index;
//    std::vector<std::map<int,int> > column_index_cache;

    /* actual compressed partitions matrices (Ai) */
    std::map<int, CompRow_Mat_double> parts;

    /* indices of the partitions owned by each masters */
    std::vector<std::vector<int> > partitionsSets;

    /* local actual compressed partitions matrices (Ai) */
    std::vector<CompRow_Mat_double> partitions;
//    std::vector<int> parts_id;

    /* number of partitions local to the master */
    int nb_local_parts;

    /**************************************************************************
     * Hybrid scheme
    **************************************************************************/
    /* The communicator shared by CG masters */
    mpi::communicator inter_comm;

    /* The communicator of local slaves */
    mpi::communicator intra_comm;

    /*! The type of process
     * - 0 the process will behave as a CG master
     * - 1 the process will be a  MUMPS Slave
     */
    int instance_type;

    /* Type of all processes */
    std::vector<int> instance_type_vect;

    /* Node of each master */
    std::vector<int> masters_node;

    /* Map containing nodes and the currently unassigned processes they contain */
    std::vector<std::vector<int>> node_map_slaves;

    /* communication map MPI-node */
    std::vector<int> mpi_map;

    /* Rank of the master for a slave */
    int my_master;

    /* Ranks of the slaves for a master */
    std::vector<int> my_slaves;

    /**************************************************************************
     * Distribution
    **************************************************************************/
    /* for each local partition, correspondance actual column-compressed column */
    std::vector<std::map<int,int> > glob_to_part;

    /* for each local partition, correspondance compressed column-actual column */
    std::vector<std::map<int,int> > part_to_glob;

    /* for the merged local partitions, correspondance actual column-compressed column */
    std::map<int,int> glob_to_local;

    /* vector of global indices from the merged local partitions */
    std::vector<int> glob_to_local_ind;

    /* defines the starting point of C in the local columns */
    int st_c_part;
    std::vector<int>::iterator st_c_part_it;
//    std::map<int,int> glob_to_local_c;

    /**
     * Contains the mutual interconnections between partitions
     * The key is the cg-master rank (in inter_comm) and the value is the column indices
     */
    std::map<int,std::vector<int> > col_interconnections;
//    std::map<int,std::vector<int> > col_inter;

    /* communication map containing 1 for process interconnected and before in the MPI ranks */
    std::vector<int> comm_map;


    /**************************************************************************
     * Cimmino
    **************************************************************************/
    /* MUMPS object to solve augmented systems */
    MUMPS mumps;

    /* size and number of non-zeros of the block diagonal augmented systems for the local partitions */
    int m_n, m_nz;
    int n_aug, nz_aug;

    /* arrays of rows/columns/values of the block diagonal augmented systems for the local partitions */
    std::vector<int> irn_aug, jcn_aug;
    std::vector<double> val_aug;

    /**************************************************************************
     * Augmentation
    **************************************************************************/
    /* indices of the kept columns in S after filtering */
    std::vector<int> selected_S_columns;

    /* indices of the removed columns in S after filtering */
    std::vector<int> skipped_S_columns;

    /* stC(i) contains the starting point of C corresponding to the first 
     * interconnected partition */
    std::vector<int> stC;

    /* MUMPS object to solve on S */
    MUMPS mumps_S;

    /* Matrix S (local) */
    Coord_Mat_double S;

/*******************************************************************************
 * Methods
 ******************************************************************************/
    /**************************************************************************
     * General Utils
    **************************************************************************/
    void initializeLog();

    /**************************************************************************
     * Scaling
    **************************************************************************/
    void scaling();
    void scaleMatrix();
    void diagScaleMatrix(std::vector<double> & , std::vector<double> & );
    void diagScaleStart (MV_ColMat_double &X, std::vector<double> &dcol);
    void diagScaleRhs (VECTOR_double &b, std::vector<double> &drow);
    void diagScaleRhs (MV_ColMat_double &B, std::vector<double> &drow);


    /**************************************************************************
     * Partitioning
    **************************************************************************/
    void partitionMatrix();
    void permutation(int *partweights_perm);
    void analyseFrame();
    void partitionWeights(std::vector<std::vector<int> > &parts,
                          std::vector<int> weights, int nb_parts);
//    void partitioning(std::vector<std::vector<int> > &, std::vector<int>, int);

    /**************************************************************************
     * Initialization
    **************************************************************************/
    void allocateMumpsSlaves(MUMPS &);
    void createAugmentedSystems(int &n_aug,
                                int &nz_aug,
                                std::vector<int> &irn_aug,
                                std::vector<int> &jcn_aug,
                                std::vector<double> &val_aug);
    void initializeDirectSolver();
    void initializeMumps(MUMPS &, bool local);
    void initializeMumps(MUMPS &);
    void analyseAugmentedSystems(MUMPS &);
    void factorizeAugmentedSystems(MUMPS &);

    /**************************************************************************
     * Distribution
    **************************************************************************/
    void createInterCommunicators();
    void distributePartitions();
    void createInterconnections();
    void distributeData();
    void distributeRhs();
//    void distributeNewRhs();

    /**************************************************************************
     * Cimmino
    **************************************************************************/
    void bcg(MV_ColMat_double &b);
    int gqr(MV_ColMat_double &P, MV_ColMat_double &AP, MV_ColMat_double &R, int s, bool use_a);
    int gqr(MV_ColMat_double &p, MV_ColMat_double &ap, MV_ColMat_double &r, CompCol_Mat_double g, int s, bool use_a);
    void gmgs2(MV_ColMat_double &p, MV_ColMat_double &ap, MV_ColMat_double &r, int s, bool use_a);
    void gmgs2(MV_ColMat_double &p, MV_ColMat_double &ap, MV_ColMat_double &r, CompCol_Mat_double g, int s, bool use_a);
    MV_ColMat_double sumProject(double alpha,
                                MV_ColMat_double &Rhs,
                                double beta,
                                MV_ColMat_double &X);

    /**************************************************************************
     * Augmentation
    **************************************************************************/
    void augmentMatrix(std::vector<CompCol_Mat_double > &loc_parts);
    void cijAugmentMatrix(std::vector<CompCol_Mat_double > &loc_parts);
    void aijAugmentMatrix(std::vector<CompCol_Mat_double > &loc_parts);
    void solveABCD(MV_ColMat_double &b);
    Coord_Mat_double buildS();
    Coord_Mat_double buildS(std::vector<int>);
    void buildS(std::vector<int> &rows,
                std::vector<int> &cols,
                std::vector<double> &vals);
    void buildS(std::vector<int> &rows,
                std::vector<int> &cols,
                std::vector<double> &vals,
                std::vector<int> &columns_to_build);
    MV_ColMat_double spSimpleProject(std::vector<int> mycols);
    void spSimpleProject(std::vector<int> mycols, std::vector<int> &vrows,
                         std::vector<int> &vcols, std::vector<double> &vvals);
    MV_ColMat_double solveS ( MV_ColMat_double &f );
    inline int S_nnz() { return S_vals.size(); }

    /* WIP */
    VECTOR_double pcgS ( VECTOR_double &b );
    MUMPS buildM();
    VECTOR_double solveM ( MUMPS &mu, VECTOR_double &z );
    MV_ColMat_double prodSv(MV_ColMat_double &);

    /**************************************************************************
     * Utils
    **************************************************************************/
    double compute_rho(MV_ColMat_double &X, MV_ColMat_double &U);
    double ddot(VECTOR_double &p, VECTOR_double &ap);
    void get_nrmres(MV_ColMat_double &x,
                    MV_ColMat_double &b,
                    VECTOR_double &nrmR,
                    VECTOR_double &nrmX);
    void centralizeVector(double *dest, int dest_lda, int dest_ncols,
                          double *src,  int src_lda,  int src_ncols,
                          std::vector<int> globalIndex, double *scale);
    void waitForSolve();
};


void configure_logger(std::string log_file);
void logger_set_filename(std::string log_file);
double or_bin(double &a, double &b);
std::vector<int> sort_indexes(const int *v, const int nb_el);

#endif // ABCD_HXX
