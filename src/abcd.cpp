#include "abcd.h"

using namespace Eigen;
using namespace std;

/// The gateway function that launches all other options
int abcd::bc(int job)
{
    mpi::communicator world;
    switch(job) {

    case -1:
        if(world.rank() == 0)
            abcd::init();
        world.barrier();
        break;

    case 1:
        if(world.rank() == 0) {
            abcd::preprocess();
            abcd::partition();
            abcd::frame_analysis();
        }
        world.barrier();
        break;

    case 2:
        // Create the group of CG instances
        abcd::inter_group_mapping();

        if(instance_type == 0) {
            abcd::distribute_parts();
        }
        break;

    case 3:
        break;

    default:
        // Wrong job id
        return -1;
    }
    return 0;
}

/// Creates and intializes the matrix object
void abcd::init()
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> elements;
    
    // Check that the matrix data is present
    if(irn == 0 ||
            jcn == 0 ||
            val == 0) {

        // Hey!! where is my data?
        throw - 1;
    }

    mtx = SparseMatrix<double>(m, n);
    if(sym) {
        for(unsigned k = 0; k < nz; ++k) {
            elements.push_back(T(irn[k] - 1, jcn[k] - 1, val[k]));
            elements.push_back(T(jcn[k] - 1, irn[k] - 1, val[k]));
        }
    } else {
        for(unsigned k = 0; k < nz; ++k) {
            elements.push_back(T(irn[k] - 1, jcn[k] - 1, val[k]));
        }
    }

    // create our object
    mtx.setFromTriplets(elements.begin(), elements.end());
    mtx.makeCompressed();
    nz = mtx.nonZeros();
}


/// Set the defaults in the constructor
abcd::abcd()
{
    start_index = 0;
}

abcd::~abcd()
{
}
