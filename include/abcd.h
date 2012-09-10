/*
 * abcd.h
 *
 *  Created on: Aug 15, 2012
 *      Author: Mohamed Zenadi
 */

#ifndef ABCD_HXX_
#define ABCD_HXX_

#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include "mmio.h"
#include "mpi.h"

#include "dmumps_c.h"

#include <Eigen/Sparse>

#include <boost/mpi.hpp>
#include <boost/range/irange.hpp>
#include <boost/foreach.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include <boost/progress.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace Eigen;
using namespace std;
using namespace boost;
using namespace boost::numeric;
using namespace ublas;


class abcd
{
private:
    // Types to be used localy

    double *absa;
    double *vnrma;

    void initialize();

    // preprocess stuffs
    /**
     * Scales the matrix
     * @norm the norm at which the matrix is scaled
     */
    void scaleMatrix(int norm);
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
     * Compresses the partitions and analyses the interconnections between them
     */
    void analyseFrame();

    // Communication stuffs
    void createInterComm();
    void distributePartitions();

    // Cimmino
    void initializeCimmino();
    
    // MUMPS
    int m_n;
    int m_nz;
    DMUMPS_STRUC_C mumps;
    void initializeMumps(bool local);
    void initializeMumps();
    void createAugmentedSystems();
    void analyseAugmentedSystems();
    void allocateMumpsSlaves();
    void factorizeAugmentedSystems();
    std::vector<int> my_slaves;
    int my_master;
    
    // MUMPS setters and getters
    inline void setMumpsIcntl(int i, int v) { mumps.icntl[ i - 1 ] = v ; }
    inline void setMumpsCntl(int i, double v) { mumps.cntl[ i - 1 ] = v ; }
    inline int getMumpsInfo(int i) { return mumps.info[ i - 1 ]; }
    inline double getMumpsRinfo(int i) { return mumps.rinfo[ i - 1 ]; }

    // SOme utilities
    void partitionWeights(std::vector<int> &, std::vector<int>, int);


    /*
     * Scaling information
     */
    VectorXd drow_;
    VectorXd dcol_;

public:
    /***************************************************************************
     * Matrix information
     */
    unsigned m;
    unsigned n;
    unsigned nz;

    /***************************************************************************
     * Temporary data about the matrix
    ***************************************************************************/
    int *irn;
    int *jcn;
    double *val;

    /***************************************************************************
     * Matrix properties
    ***************************************************************************/
    bool sym; /// Symmetry
    int start_index; /// To define wether it's Fortran-Style (1) or C-Style (0)

    /***************************************************************************
     * The matrix object itself
    ***************************************************************************/
    Eigen::SparseMatrix<double, RowMajor> mtx;
    std::vector<Eigen::SparseMatrix<double, RowMajor> > parts;

    /***************************************************************************
     * Partitioning informations
    ***************************************************************************/

    /**
     * Defines the type of partitioning :
     * * 1 for manual partitioning (nbparts and nbrows are suplied)
     * * 2 for automatic partitioning (nbparts only is needed)
     */
    int partitioning_type;
    int nbparts; /// The number of partitions
    ArrayXi strow; /// The starting row index of each partition
    ArrayXi nbrows; /// The number of rows per partition
    /// A reverse index of columns, contains the original index of each column for each partition
    std::vector<std::vector<int> > columns_index;
    /// A merge of col_index vectors, determines non-null columns in all local partitions
    std::vector<int> global_column_index;
    /**
     * Contains the mutual interconnections between partitions
     * The key is the cg-master rank (in inter_comm) and the value is the column indices
     */
    std::map<int,std::vector<int> > col_interconnections;
    /// Contains the partitions that are handled by this instance
    std::vector<int> parts_id;



    /***************************************************************************
     * Communication info
    ***************************************************************************/
    int instance_type; /// The type of process : 0 for CG, 1 for MUMPS Slave
    int parallel_cg; /// The number of parallel CG instances

    mpi::communicator inter_comm; /// The communicator shared by CG masters
    mpi::communicator intra_comm; /// The communicator of local slaves


    int icntl[10];
    double dcntl[10];


    int bc(int);
    void preprocess();
    abcd();
    ~abcd();

};

typedef std::pair<double,int> dipair;
bool ip_comp(const dipair &, const dipair &);

#endif // ABCD_HXX
