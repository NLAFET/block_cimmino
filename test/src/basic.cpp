#include <abcd.h>

using namespace boost;
int main()
{
    mpi::communicator world;

    abcd obj;

    if (world.rank() == 0) {
        obj.sym = true;
        obj.n = 5;
        obj.m = 5;
        obj.nz = 10;

        obj.irn = new int[12];
        obj.jcn = new int[12];
        obj.val = new double[12];
    } else {
    }
    return 0;
}
