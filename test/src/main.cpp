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

        /* read the file and its content */
        //f = fopen("/home/knuthy/work/stash/offshore/offshore_pr.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/gre_1107/gre_1107_pr.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/gre_1107/gre_1107.mtx", "r");
        f = fopen("/home/knuthy/work/stash/bayer01/bayer01_pr.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/b1_ss/b1_ss.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/SiO2/SiO2_pr.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/SiO/SiO.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/LFAT5/LFAT5.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/bone010/bone010_pr.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/bayer01/bayer01.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/RM07R/RM07R.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/GT01R/GT01R.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/lhr34c/lhr34c.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/bayer01/bayer01.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/ASIC_320ks/ASIC_320ks.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/PR02R/PR02R.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/pores_3/pores_3.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/s3_sym.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/Huhs_DNA_CC.mtx", "r");
        //f = fopen("/home/knuthy/work/stash/1970", "r");
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
        obj.icntl[10]= 1;

        obj.nbparts = 16;

        obj.partitioning_type = 1;
        obj.nbrows = VECTOR_int(obj.nbparts);
        //gre_1107
        //obj.nbrows << 219, 218, 223, 224, 223;
        //
        //bayer01
        int nr[] = {3659, 3653, 3610, 3662, 3582, 3583, 3634, 3623, 3604, 3569, 3608, 3588, 3589, 3603, 3569, 3599};
        obj.nbrows = VECTOR_int(nr, obj.nbparts);

        /*sio2*/
        //obj.nbrows << 5117,5213,5200,5158,5192,5205,5145,5137,5136,5181,5181,5186,5245,5193,5249,5112,5203,5195,5141,5198,5140,5200,5146,5145,5146,5145,5212,5151,5206,5253;

        /* offshore */
        //obj.nbrows <<   8613,8761,8761,8669,8754,8710,8710,8669,8668,8732,8775,8734, 8735,8687,8686,8666,8656,8656,8625,8624,8694,8694,8519,8519, 8609,8608,8564,8564,8564,8563;
        //obj.nbrows << 10967,11035,10957,10907,10908,11038,11079,10993,11083,11014,11013,11022,11010,11016,11004,11004,10915,11027,11028,10993,10996,10980,10934,10934,11011,10941,10942,11018,11067,11161,11017,11017,11016,10899,10994,10992,11069,10982,11019,11013,11028,11031,10972,10956,10955,10869,10869,10964,10961,10961,11061,11035,10948,10960,11042,11042,10986,10987,10967,10965,10966,10905,10982,11034,10907,11007,11007,10865,10900,10902,10862,10880,10832,10781,10781,10842,10776,10776,10958,10911,10912,10873,10961,10947,10877,10952,10952,10926,10927,10928;

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
            obj.itmax = 10;
            // works only in sequential for the moment
            obj.use_xf = true;
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
