#include "abcd.h"
#include <iostream>
#include <cstdio>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Test_abcd Unit Tests"

#include <boost/test/unit_test.hpp>

using namespace std;

BOOST_AUTO_TEST_CASE(my_test)
{
    BOOST_REQUIRE(true);
    
    mpi::environment env();
    mpi::communicator world;
    
    cout << "Hey, I'm process " << world.rank() << " of " << world.size() << endl;
    
    abcd obj;
    FILE *f;
    MM_typecode mat_code;

    // read the file and its content
    f = fopen("/home/knuthy/stash/gre_1107/gre_1107.mtx", "r");
//     f = fopen("/home/knuthy/stash/s3_sym.mtx", "r");
//     f = fopen("/home/knuthy/stash/Huhs_DNA_CC.mtx", "r");
//     f = fopen("/home/knuthy/stash/1970", "r");
    mm_read_banner(f, &mat_code);
    mm_read_mtx_crd_size(f, (int *)&obj.m, (int *)&obj.n, (int *)&obj.nz);

    if (mm_is_symmetric(mat_code))
        obj.sym = true;
    else
        obj.sym = false;

    // allocate the arrays
    obj.irn = new int[obj.nz];
    obj.jcn = new int[obj.nz];
    obj.val = new double[obj.nz];

    obj.start_index = 1;
    obj.icntl[9] = 2;
    
    obj.nbparts = 5;
    obj.partitioning_type = 2;

    mm_read_mtx_crd_data(f, obj.m, obj.n, obj.nz, obj.irn, obj.jcn, obj.val,
                         mat_code);

    cout << "Matrix information : ";
    cout << "m = " << obj.m << "; n = " << obj.n << "; nz = " << obj.nz << endl;

    fclose(f);

    BOOST_CHECK(obj.bc(-1) == 0);
    BOOST_CHECK(obj.bc(1) == 0);

}