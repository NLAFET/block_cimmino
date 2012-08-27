#include <iostream>
#include <string>

#include <cstdio>

#include "abcd.h"
#include "mpi.h"


using namespace std;

int main(int argc, char *argv[])
{
    abcd obj;
    FILE *f;
    MM_typecode mat_code;
    
    MPI_Init();

    // read the file and its content
//     f = fopen("/home/knuthy/stash/gre_1107/gre_1107.mtx", "r");
    f = fopen("/home/knuthy/stash/s3_sym.mtx", "r");
    mm_read_banner(f, &mat_code);
    mm_read_mtx_crd_size(f, &obj.m, &obj.n, &obj.nz);

    if (mm_is_symmetric(mat_code))
        obj.sym = true;
    else
        obj.sym = false;

    // allocate the arrays
    obj.irn = new int[obj.nz];
    obj.jcn = new int[obj.nz];
    obj.val = new double[obj.nz];

    obj.start_index = 1;

    mm_read_mtx_crd_data(f, obj.m, obj.n, obj.nz, obj.irn, obj.jcn, obj.val,
                         mat_code);

    cout << "Matrix information : ";
    cout << "m = " << obj.m << "; n = " << obj.n << "; nz = " << obj.nz << endl;

    fclose(f);

    obj.icntl[9] = 2;
    obj.bc(-1);

    obj.bc(1);

    return 0;
}
