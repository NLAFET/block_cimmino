#include <abcd.h>

using namespace Eigen;
using namespace std;

/// The gateway function that launches all other options
int abcd::bc(int job)
{
    mpi::communicator world;

    switch(job) {

    case -1:
        if(world.rank() == 0) {
            abcd::initialize();
        }
        world.barrier();
        break;

    case 1:
        if(world.rank() == 0) {
            abcd::preprocess();
            abcd::partitionMatrix();
            abcd::analyseFrame();
        }
        world.barrier();
        break;

    case 2:
        // Create the group of CG instances
        abcd::createInterComm();

        if(instance_type == 0) {
            abcd::distributePartitions();
        }
        abcd::initializeCimmino();
        if(instance_type == 0)
            cout << "[+] Launching MUMPS factorization" << endl;
        abcd::factorizeAugmentedSystems();
        break;

    case 3:
        if(instance_type == 0) {
            inter_comm.barrier();
            abcd::distributeRhs();
            abcd::bcg();
        }
        world.barrier();
        break;

    default:
        // Wrong job id
        return -1;
    }
    return 0;
}

/// Creates and intializes the matrix object
void abcd::initialize()
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
            if(irn[k] != jcn[k])
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
    nrhs = 1;
    start_index = 0;
    block_size = 1;
    use_xk = false;
}

abcd::~abcd()
{
}
