#include <iostream>
#include <string>
#include <cstdio>

#include "abcd.h"

using namespace std;

int main(int argc, char *argv[])
{
    mpi::environment env;
    mpi::communicator world;

    // This should be done only by the master
    if(world.rank() == 0) {
        abcd obj;
        FILE *f;
        MM_typecode mat_code;

        // read the file and its content
        //f = fopen("/home/knuthy/stash/gre_1107/gre_1107.mtx", "r");
//         f = fopen("/home/knuthy/stash/lhr34c/lhr34c.mtx", "r");
//        f = fopen("/home/knuthy/stash/bayer01/bayer01.mtx", "r");
//         f = fopen("/home/knuthy/stash/ASIC_320ks/ASIC_320ks.mtx", "r");
//         f = fopen("/home/knuthy/stash/PR02R/PR02R.mtx", "r");
//         f = fopen("/home/knuthy/stash/pores_3/pores_3.mtx", "r");
//         f = fopen("/home/knuthy/stash/s3_sym.mtx", "r");
//         f = fopen("/home/knuthy/stash/Huhs_DNA_CC.mtx", "r");
         f = fopen("/home/knuthy/stash/1970", "r");
        mm_read_banner(f, &mat_code);
        mm_read_mtx_crd_size(f, (int *)&obj.m, (int *)&obj.n, (int *)&obj.nz);

        if(mm_is_symmetric(mat_code))
            obj.sym = true;
        else
            obj.sym = false;

        // allocate the arrays
        obj.irn = new int[obj.nz];
        obj.jcn = new int[obj.nz];
        obj.val = new double[obj.nz];

        obj.start_index = 1;
        obj.icntl[9] = 0;

        obj.nbparts = 4;
        obj.partitioning_type = 2;

        obj.parallel_cg = obj.nbparts < world.size() ? obj.nbparts : world.size();

        mm_read_mtx_crd_data(f, obj.m, obj.n, obj.nz, obj.irn, obj.jcn, obj.val,
                             mat_code);

        cout << "Matrix information : ";
        cout << "m = " << obj.m << "; n = " << obj.n << "; nz = " << obj.nz << endl;

        fclose(f);

        try {
            obj.bc(-1);
            obj.bc(1);
            obj.bc(2);

            obj.block_size = 2;
            obj.itmax = 1000;
            obj.rhs = new double[obj.m_l * obj.nrhs];
            for(int j = 0; j < obj.nrhs; j++)
                for(int i = 0; i < obj.m_l; i++)
                    obj.rhs[i + j * obj.m_l] = j+1;

            obj.bc(3);
        } catch(int e) {
            cout << "Error code : " << e << endl;
        }
    } else {
        abcd obj;

        try {
            obj.bc(-1);
            obj.bc(1);
            obj.bc(2);
            obj.bc(3);
        } catch(int e) {
            cout << "Error code : " << e << endl;
        }
    }
    world.barrier();


    return 0;
}
