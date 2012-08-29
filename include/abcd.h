/*
 * abcd.h
 *
 *  Created on: Aug 15, 2012
 *      Author: Mohamed Zenadi
 */

#ifndef ABCD_HXX_
#define ABCD_HXX_

#include <iostream>
#include <vector>
#include <cstdio>
#include "mmio.h"
#include "mpi.h"

#include <Eigen/Sparse>

#include <boost/mpi.hpp>
#include <boost/range/irange.hpp>
#include <boost/foreach.hpp>
#include <boost/range/algorithm.hpp>

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

    void init();

    // preprocess stuffs
    /**
     * Scales the matrix
     * @norm the norm at which the matrix is scaled
     */
    void scal_matrix(int norm);
    /**
     * Computes the norm of the matrix
     * @todo implement it!
     */
    void comp_norm();

    // structure functions
    /// Partitions the matrix into abcd::nbrows
    void partition();
    /**
     * Analyses the structure of each partition
     * Compresses the partitions and analyses the interconnections between them
     */
    void frame_analysis();

    // Communication stuffs
    void inter_group_mapping();
    void distribute_parts();


    // SOme utilities
    void partition_parts(std::vector<std::vector<int> > &, std::vector<int>);


    /*
     * Scaling information
     */
    VectorXd drow;
    VectorXd dcol;

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
    std::vector<Eigen::SparseMatrix<double, ColMajor> > parts;

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
    std::vector<std::vector<int> > col_index;
    /// Contains the mutual interconnections between partitions
    std::vector<std::vector<int> > col_interconnections;



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

#endif // ABCD_HXX
