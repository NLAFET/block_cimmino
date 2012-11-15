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
        //f = fopen("/home/knuthy/work/stash/gre_1107/gre_1107_pr.mtx", "r");
         //f = fopen("/home/knuthy/work/stash/gre_1107/gre_1107.mtx", "r");
        f = fopen("/home/knuthy/work/stash/bayer01/bayer01_pr.mtx", "r");
       //f = fopen("/home/knuthy/work/stash/bayer01/bayer01.mtx", "r");
//        f = fopen("/home/knuthy/work/stash/RM07R/RM07R.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/GT01R/GT01R.mtx", "r");
//        f = fopen("/home/knuthy/work/stash/lhr34c/lhr34c.mtx", "r");
//        f = fopen("/home/knuthy/work/stash/bayer01/bayer01.mtx", "r");
//        f = fopen("/home/knuthy/work/stash/ASIC_320ks/ASIC_320ks.mtx", "r");
//        f = fopen("/home/knuthy/work/stash/PR02R/PR02R.mtx", "r");
//        f = fopen("/home/knuthy/work/stash/pores_3/pores_3.mtx", "r");
//        f = fopen("/home/knuthy/work/stash/s3_sym.mtx", "r");
//        f = fopen("/home/knuthy/work/stash/Huhs_DNA_CC.mtx", "r");
//        f = fopen("/home/knuthy/work/stash/1970", "r");
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
        obj.icntl[9] = 2;

        obj.nbparts = 16;

        obj.partitioning_type = 1;
        obj.nbrows = ArrayXi(obj.nbparts);
        //obj.nbrows << 219, 218, 223, 224, 223;
        obj.nbrows << 3659, 3653, 3610, 3662, 3582, 3583, 3634, 3623, 3604, 3569, 3608, 3588, 3589, 3603, 3569, 3599;

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

            obj.block_size = 1;
            obj.itmax = 1000;
            obj.use_xf = false;
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
