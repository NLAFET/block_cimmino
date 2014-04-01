// the solver's header file
#include "abcd.h"

// use boost::mpi for simplicity, the user can use which ever he wants
#include "mpi.h"
#include <boost/mpi.hpp>

int main(int argc, char* argv[]) 
{
    // equivalent to MPI_Initialize
    mpi::environment env(argc, argv);

    // obtain the WORLD communicator, by default the solver uses it
    mpi::communicator world;

    // create one instance of the abcd solver per mpi-process
    abcd obj;

    if(world.rank() == 0) { // the master
        // we create a 5x5 matrix for a 1D mesh + three-point stencil 
        obj.m = 10; // number of rows
        obj.n = obj.m; // number of columns
        obj.nz = 2*obj.m - 1; // number of nnz in the lower-triangular part
        obj.sym = true;

        // allocate the arrays
        obj.irn = new int[obj.nz];
        obj.jcn = new int[obj.nz];
        obj.val = new double[obj.nz];

        // initialize the matrix
        // Note: the matrix is stored in 1-based format
        size_t pos = 0;
        for (size_t i = 1; i < obj.m; i++) {

            // the diagonal
            obj.irn[pos] = i;
            obj.jcn[pos] = i;
            obj.val[pos] = 2.0;

            pos++;

            // the lower-triangular part
            obj.irn[pos] = i + 1;
            obj.jcn[pos] = i;
            obj.val[pos] = -1.0;

            pos++;
        }

        // the last diagonal element
        obj.irn[pos] = obj.m;
        obj.jcn[pos] = obj.m;
        obj.val[pos] = 2.0;

        pos++;

        // set the rhs
        obj.rhs = new double[obj.m];
        for (size_t i = 0; i < obj.m; i++) {
            obj.rhs[i] = ((double) i + 1)/obj.m;
        }

        // this should be correct
        assert(pos == obj.nz);

        // ask the solver to guess the number of partitions
        obj.icntl[Controls::part_guess] = 1;

        
    } else { // the workers

    }
    try {
        obj(-1);
        obj(5); // equivalent to running 1, 2 and 3 successively

        if (world.rank() == 0) {
            cout << obj.sol << endl;
        }
    } catch (runtime_error err) {
        cout << "An error occured: " << err.what() << endl;
    }

  return 0;
}
