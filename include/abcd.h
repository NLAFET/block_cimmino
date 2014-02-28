/*
 * abcd.h
 *
 *  Created on: Aug 15, 2012
 *      Author: Mohamed Zenadi
 */

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
#include "mmio.h"

#include <mumps.h>
#include "dmumps_c.h"

#include <boost/mpi.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/progress.hpp>
#include <boost/lambda/lambda.hpp>
//#include <boost/range/adaptors.hpp>
//#include <boost/range/algorithm.hpp>

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

/* Some macros*/
#define IRANK inter_comm.rank()
#define IBARRIER inter_comm.barrier()

#define IFMASTER if(inter_comm.rank() == 0)
#define TIC t = MPI_Wtime()
#define TOC MPI_Wtime() - t

#define inter_rank() inter_comm.rank()
#define inter_barrier() inter_comm.barrier()
#define inter_master() if(inter_comm.rank() == 0)

using namespace std;
using namespace boost;
using namespace boost::numeric;
//using namespace boost::adaptors;

class abcd
{
private:
    // Types to be used localy
    double nrmA;
    double nrmB;
    double nrmXf;
    double nrmMtx;

    // preprocess stuffs
    void preprocess();
    /**
     * Scales the matrix
     * @norm the norm at which the matrix is scaled
     */
    void scaleMatrix(int norm);
    void diagScaleMatrix(VECTOR_double , VECTOR_double );
    void diagScaleRhs(VECTOR_double &);
    void diagScaleRhs(MV_ColMat_double &);
    /**
     * Computes the norm of the matrix
     * @todo implement it!
     */
    void computeNorms();

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
    
    /*-----------------------------------------------------------------------------
     *  Build the augmented version of the matrix (ABCD)
     *-----------------------------------------------------------------------------*/
    void augmentMatrix(std::vector<CompCol_Mat_double > &loc_parts);

    // Communication stuffs
    void createInterComm();
    void distributePartitions();
    void createInterconnections();

    void distributeData();

    // Cimmino
    void initializeCimmino();
    void distributeRhs();
    void distributeNewRhs();
    void bcg(MV_ColMat_double &b);
    void solveABCD(MV_ColMat_double &b);
    MV_ColMat_double solveS ( MV_ColMat_double &f );
    Coord_Mat_double buildS();
    Coord_Mat_double buildS(std::vector<int>);
    MUMPS buildM();
    VECTOR_double solveM ( MUMPS &mu, VECTOR_double &z );
    MV_ColMat_double prodSv(MV_ColMat_double &);
    VECTOR_double pcgS ( VECTOR_double &b );
    std::vector<int> selected_S_columns;
    std::vector<int> skipped_S_columns;
    int gqr(MV_ColMat_double &P, MV_ColMat_double &AP, MV_ColMat_double &R, int s, bool use_a);
    int gqr(MV_ColMat_double &p, MV_ColMat_double &ap, MV_ColMat_double &r, CompCol_Mat_double g, int s, bool use_a);
    void gmgs(MV_ColMat_double &p, MV_ColMat_double &ap, MV_ColMat_double &r, int s, bool use_a);
    void gmgs(MV_ColMat_double &p, MV_ColMat_double &ap, MV_ColMat_double &r, CompCol_Mat_double g, int s, bool use_a);
    double compute_rho(MV_ColMat_double &X, MV_ColMat_double &U, double thresh);
    std::vector<double> normres;
    int size_c;

    // MUMPS
    int m_n;
    int m_nz;
    MUMPS mumps;
    void initializeMumps(MUMPS &, bool local);
    void initializeMumps(MUMPS &);
    void createAugmentedSystems(MUMPS &);
    void analyseAugmentedSystems(MUMPS &);
    void allocateMumpsSlaves(MUMPS &);
    void factorizeAugmentedSystems(MUMPS &);
    std::vector<int> my_slaves;
    int my_master;
    MV_ColMat_double sumProject(double alpha, MV_ColMat_double &Rhs, double beta, MV_ColMat_double &X);
    MV_ColMat_double coupleSumProject(double alpha, MV_ColMat_double &Rhs, double beta, MV_ColMat_double &X, int my_bro);

    MV_ColMat_double simpleProject(MV_ColMat_double &X);
    MV_ColMat_double spSimpleProject(std::vector<int> mycols);
    void spSimpleProject(std::vector<int> mycols, std::vector<int> &vrows,
            std::vector<int> &vcols, std::vector<double> &vvals);

    void waitForSolve();
    std::vector<int> comm_map;

    // SOme utilities
    void partitionWeights(std::vector<std::vector<int> > &, std::vector<int>, int);
    void partitioning(std::vector<std::vector<int> > &, std::vector<int>, int);
    double ddot(VECTOR_double &p, VECTOR_double &ap);
    void get_nrmres(MV_ColMat_double &x, MV_ColMat_double &b, double &nrmR, double &nrmX, double &nrmXfmX);


    /*
     * Scaling information
     */
    VECTOR_double drow_;
    VECTOR_double dcol_;

    /***************************************************************************
     * The matrix object itself
    ***************************************************************************/
    Coord_Mat_double S;
    std::map<int, CompRow_Mat_double> parts;
    std::vector<CompRow_Mat_double> partitions;
    std::vector<vector<int> > p_sets;

    MV_ColMat_double Xf;
    MV_ColMat_double B;
    MV_ColMat_double RRand;

    CompRow_Mat_double A;
    std::vector<int> row_perm;

    bool runSolveS;

public:
    /***************************************************************************
     * Matrix information
     */
    int m;
    int n;
    int nz;
    int nrhs;
    int n_l, m_l, nz_l;
    int n_o, m_o, nz_o;

    //! The path where to write the matrix \f$PD_rAD_cP^T\f$
    std::string write_problem;
    //! The path where to write the matrix \f$S\f$ 
    std::string write_s;


    /***************************************************************************
     * Temporary data about the matrix
    ***************************************************************************/
    int *irn;
    int *jcn;
    double *val;
    double *rhs;

    MV_ColMat_double Xk;
    MV_ColMat_double sol;
    

    /***************************************************************************
     * Matrix properties
    ***************************************************************************/
    bool sym; /// Symmetry
    int start_index; /// To define wether it's Fortran-Style (1) or C-Style (0)

    /***************************************************************************
     * Partitioning informations
    ***************************************************************************/

    /*!
     * Defines the type of partitioning :
     * * 1 for manual partitioning (nbparts and nbrows are suplied)
     * * 2 for automatic partitioning (nbparts only is needed)
     */
    int partitioning_type;

    /*!
     * Guess the number of partitions
     * * 1 based on the number of rows
     */
    int guessPartitionsNumber;

    int nbparts; /// The number of partitions
    int nb_local_parts;
    VECTOR_int strow; /// The starting row index of each partition
    VECTOR_int nbrows; /// The number of rows per partition
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

    /**
     * Contains the mutual interconnections between partitions
     * The key is the cg-master rank (in inter_comm) and the value is the column indices
     */
    std::map<int,std::vector<int> > col_interconnections;
    std::map<int,std::vector<int> > col_inter;
    /// Contains the partitions that are handled by this instance
    std::vector<int> parts_id;

    int block_size;
    int itmax;
    double threshold;

    bool verbose;


    /***************************************************************************
     * Communication info
    ***************************************************************************/
    /// The number of parallel CG instances
    int parallel_cg;

    /// The communicator shared by CG masters
    mpi::communicator inter_comm;
    /// The communicator of local slaves
    mpi::communicator intra_comm; 


    /*! Control parameters
     *
     *  - __icntl[10]__: augmentation type
     *      - 0 : No augmentation, classical Block Cimmino will be used
     *      - 1 : \f$C_{i,j}/-I\f$ based augmentation
     *      - 2 : \f$A_{i,j}/-A_{j,i}\f$ based augmentation
     *  - __icntl[11]__: augmentation analysis
     *  - __icntl[12]__: compute \f$\bar{A}^+w\f$ only
     *  - __icntl[13]__: use dense right-hand sides when building \f$S\f$
     *  - __icntl[14]__: the blocking-factor used when building \f$S\f$
     *  - __icntl[15]__: solve \f$Sz =f\f$ iteratively
     */
    int icntl[20];
    double dcntl[20];
    int info[20];
    double dinfo[20];

    int initializeMatrix();
    int preprocessMatrix();
    int factorizeAugmentedSystems();
    int solveSystem();

    int bc(int);
    abcd();
    ~abcd();

};

typedef std::pair<double,int> dipair;
bool ip_comp(const dipair &, const dipair &);
template <class K, class V> std::vector<K> get_keys(std::map<K,V> my_map);
double or_bin(double &a, double &b);
void setVal(int *lst, int sz, int ival);
vector<int> sort_indexes(const int *v, const int nb_el);
template <typename T> vector<int> sort_indexes(const vector<T> &v);

#endif // ABCD_HXX
