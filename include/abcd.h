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

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

using namespace Eigen;
using namespace std;
namespace mpi = boost::mpi;

class abcd
{
private:
    // Types to be used localy
    
    double *absa;
    double *vnrma;

    int init();

    // preprocess stuffs
    int scal_matrix(int);
    int comp_norm();
    
    // structure functions
    int partition();
    
    // Communication stuffs
    int inter_group_mapping();
    
    /*!
     * Scaling information
     */
    VectorXd drow;
    VectorXd dcol;

public:
    /*
     * Matrix information
     */
    unsigned m;
    unsigned n;
    unsigned nz;

    /*
     * Temporary data about the matrix
     */
    int *irn;
    int *jcn;
    double *val;

    /*
     * Matrix properties
     */
    bool sym; /// Symmetry
    int start_index; /// To define wether it's Fortran-Style (1) or C-Style (0)

    /*
     * The matrix object itself
     */
    Eigen::SparseMatrix<double> mtx;
    
    /*
     * Partitioning informations
     */
    /*!
     * Defines the type of partitioning :
     * * 1 for manual partitioning (nbparts and nbrows are suplied)
     * * 2 for automatic partitioning (nbparts only is needed)
     */
    int partitioning_type; 
    int nbparts; /// The number of partitions
    ArrayXi strow; /// The starting row index of each partition
    ArrayXi nbrows; /// The number of rows per partition


    int icntl[10];
    double dcntl[10];


    int bc(int);
    int preprocess();
    abcd();
    ~abcd();

};

#endif // ABCD_HXX
