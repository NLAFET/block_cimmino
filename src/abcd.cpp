#include "abcd.h"

using namespace Eigen;
using namespace std;

/// Creates the matrix object and initialize default parameters
int abcd::init()
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> elements;

    // Check that the matrix data is present
    if (irn == 0 ||
            jcn == 0 ||
            val == 0) {

        // Hey!! where is my data?
        return -1;
    }

    mtx = SparseMatrix<double>(m, n);
    if (sym) {
        for (unsigned k = 0; k < nz; ++k) {
            elements.push_back(T(irn[k] - 1, jcn[k] - 1, val[k]));
            elements.push_back(T(jcn[k] - 1, irn[k] - 1, val[k]));

//             mtx(jcn[k] - 1 , irn[k] - 1) = val[k];
        }
    } else {
        for (unsigned k = 0; k < nz; ++k) {
            elements.push_back(T(irn[k] - 1, jcn[k] - 1, val[k]));
        }
    }

    // create our object
    mtx.setFromTriplets(elements.begin(), elements.end());
    nz = mtx.nonZeros();

    return 0;
}

/// The gateway function that launches all other options
int abcd::bc(int job)
{
    int ret = 0;
    switch (job) {

    case -1:
        ret = abcd::init();
        if (ret != 0) return ret;
        break;

    case 1:
        ret = abcd::preprocess();
        if (ret != 0) return ret;
        ret = abcd::partition();
        if (ret != 0) return ret;
        break;

    case 2:
        ret = abcd::inter_group_mapping();
        break;

    case 3:
        break;

    default:
        // Wrong job id
        return -1;
    }

    // everything is alright, just return 0
    return ret;
}

/// Set the defaults in the constructor
abcd::abcd()
{
    start_index = 0;
}

/// The simplest destructor ever
abcd::~abcd()
{
    // Deallocate matrix data

    delete irn;
    delete jcn;
    delete val;
}
