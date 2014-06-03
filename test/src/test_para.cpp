#include "abcd.h"
#include "defaults.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "vect_utils.h"
#include <stdlib.h>

using ::testing::AtLeast;
using ::testing::Return;
using ::testing::Eq;
using ::testing::Ne;
using ::testing::Le;
using ::testing::Lt;
using namespace Controls;

#include "mpi.h"
#include <boost/mpi.hpp>

// A simple matrix generator for a regular 2D mesh + 5-point stencil 
void init_2d_lap(int m, int n, int nz, int *irn, int *jcn, double *val, int mesh_size);
void init_2d_lap(abcd &o, int mesh_size);


class AbcdTest : public ::testing::Test {
protected:
  mpi::communicator world;
  abcd obj;
};

TEST_F (AbcdTest, defaults)
{
  // default icntl
  int ic[] = {0,4,2,0,1,2,1000,1,0,0,0,256,0,0,0,0,0,0,0,0};
  std::vector<int> default_icntl(ic, ic + 20);

  EXPECT_THAT(obj.icntl, Eq(default_icntl));
}

TEST_F (AbcdTest, WrongJobOrder) 
{
  EXPECT_THROW(obj(1), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-2));

  EXPECT_THROW(obj(2), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-2));

  EXPECT_THROW(obj(3), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-2));

  EXPECT_THROW(obj(4), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-2));

  EXPECT_THROW(obj(5), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-2));

  EXPECT_THROW(obj(6), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-2));

  EXPECT_THROW(obj(8), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-1));
}

TEST_F (AbcdTest, MatrixInit) 
{
  int mesh_size = 10;
  init_2d_lap(obj, mesh_size);

  if (world.rank() == 0) {
    // number of non-zeros before expansion of the matrix
    EXPECT_THAT(obj.nz, Eq(3 * mesh_size * mesh_size - 2 * mesh_size));
  }

  ASSERT_NO_THROW(obj(-1));
  EXPECT_THAT(obj.info[Controls::status], Eq(0));

  // rerun
  EXPECT_THROW(obj(-1), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-2));

  if (world.rank() == 0) {
    // number of non-zeros after expansion of the matrix
    EXPECT_THAT(obj.nz, Eq(5 * mesh_size * mesh_size - 4 * mesh_size));
  }

  obj.icntl[Controls::part_guess] = 0;
  obj.icntl[Controls::part_type] = 1;
  EXPECT_ANY_THROW(obj(1));
  if (world.rank() == 0) {
    EXPECT_THAT(obj.info[Controls::status], Eq(-4));
  }

  obj.icntl[Controls::nbparts] = 0;
  obj.icntl[Controls::part_type] = 2;
  EXPECT_ANY_THROW(obj(1));
  if (world.rank() == 0) {
    EXPECT_THAT(obj.info[Controls::status], Eq(-3));
  }

  // more than what It should
  obj.icntl[Controls::nbparts] = 101;
  EXPECT_NO_THROW(obj(1));
  if (world.rank() == 0) {
    EXPECT_THAT(obj.icntl[Controls::nbparts], Eq(100));
  }
}

TEST_F (AbcdTest, OneSystemOneBased)
{
  obj.m = 1;
  obj.n = 1;
  obj.nz = 1;

  obj.irn = new int[1];
  obj.jcn = new int[1];
  obj.val = new double[1];
  obj.rhs = new double[1];

  obj.irn[0] = 1;
  obj.jcn[0] = 1;
  obj.val[0] = 1;
  obj.rhs[0] = 2;
  obj.nrhs = 1;
  obj.sym = false;

  EXPECT_NO_THROW(obj(-1));
  EXPECT_NO_THROW(obj(1));
  EXPECT_ANY_THROW(obj(1));
  EXPECT_NO_THROW(obj(5));

  if (world.rank() == 0) {
    EXPECT_THAT(obj.sol[0], Eq(obj.rhs[0]));
  }

  // auto generate a new RHS
  if (world.rank() == 0) {
    delete[] obj.rhs;
    obj.rhs = nullptr;
  }

  EXPECT_NO_THROW(obj(3));
  if (world.rank() == 0) {
    EXPECT_THAT(obj.sol[0], Eq(obj.rhs[0]));
  }
    

}

TEST_F (AbcdTest, OneSystemZeroBased)
{
  obj.m = 1;
  obj.n = 1;
  obj.nz = 1;

  obj.irn = new int[1];
  obj.jcn = new int[1];
  obj.val = new double[1];
  obj.rhs = new double[1];

  obj.irn[0] = 0;
  obj.jcn[0] = 0;
  obj.val[0] = 1;
  obj.rhs[0] = 2;
  obj.nrhs = 1;

  // write the problem
  obj.write_problem = "/tmp/test_file_abcd";

  EXPECT_NO_THROW(obj(-1));
  EXPECT_NO_THROW(obj(1));
  EXPECT_ANY_THROW(obj(1));
  EXPECT_NO_THROW(obj(5));

  if (world.rank() == 0) {
    EXPECT_THAT(obj.sol[0], Eq(obj.rhs[0]));
  }
}

TEST_F (AbcdTest, ClassicalCG) 
{
  init_2d_lap(obj, 10);

  EXPECT_NO_THROW( obj(-1) );

  EXPECT_NO_THROW( obj(1) );

  // go to past
  EXPECT_THROW(obj(-1), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-2));

  EXPECT_NO_THROW( obj(5) );
  if (world.rank() == 0) EXPECT_THAT(obj.info[Controls::nb_iter], Le(obj.icntl[Controls::itmax]));

  // go to wrong past
  EXPECT_THROW(obj(4), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-2));
  // go to wrong past
  EXPECT_THROW(obj(2), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-2));

  // re-launch solve
  EXPECT_NO_THROW( obj(3) );

  // no more than what we've requested
  obj.icntl[Controls::itmax] = 1;
  EXPECT_NO_THROW( obj(3) );
  if (world.rank() == 0) EXPECT_THAT(obj.info[Controls::nb_iter], Eq(obj.icntl[Controls::itmax]));

  // no negative iteration counts
  obj.icntl[Controls::itmax] = -1;
  EXPECT_THROW(obj(3), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-11));
  obj.icntl[Controls::itmax] = 1000;

  // block-size > 0
  obj.icntl[Controls::block_size] = 0;
  EXPECT_THROW(obj(3), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-10));
  obj.icntl[Controls::block_size] = -1;
  EXPECT_THROW(obj(3), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-10));
  obj.icntl[Controls::block_size] = 1;

  // rhs number > 0
  obj.nrhs = 0;
  EXPECT_THROW(obj(3), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-9));
  obj.nrhs = -1;
  EXPECT_THROW(obj(3), runtime_error);
  EXPECT_THAT(obj.info[Controls::status], Eq(-9));
  obj.nrhs = 1;
}

TEST_F (AbcdTest, BlockCG) 
{
  init_2d_lap(obj, 100);
  obj.icntl[block_size] = 4;

  EXPECT_NO_THROW(obj(-1));
  EXPECT_NO_THROW(obj(6));
}

TEST_F (AbcdTest, BlockCG_GMGS2) 
{
  init_2d_lap(obj, 100);
  obj.icntl[block_size] = 4;
  obj.icntl[use_gmgs2] = 1;

  EXPECT_NO_THROW(obj(-1));
  EXPECT_NO_THROW(obj(6));
}

TEST_F (AbcdTest, BlockCG_PaToH_0) 
{
  init_2d_lap(obj, 100);
  obj.icntl[block_size] = 4;
  obj.icntl[part_type] = 3;
  obj.dcntl[part_imbalance] = 0;

  EXPECT_NO_THROW(obj(-1));
  EXPECT_NO_THROW(obj(6));
}

TEST_F (AbcdTest, BlockCG_PaToH_05) 
{
  init_2d_lap(obj, 100);
  obj.icntl[block_size] = 4;
  obj.icntl[part_type] = 3;
  obj.dcntl[part_imbalance] = 0.5;

  EXPECT_NO_THROW(obj(-1));
  EXPECT_NO_THROW(obj(6));
}

TEST_F (AbcdTest, BlockCG_PaToH_15) 
{
  init_2d_lap(obj, 100);
  obj.icntl[block_size] = 4;
  obj.icntl[part_type] = 3;
  obj.dcntl[part_imbalance] = 1.5;

  EXPECT_NO_THROW(obj(-1));
  EXPECT_NO_THROW(obj(6));
}

TEST_F (AbcdTest, BlockCG_PaToH_25) 
{
  init_2d_lap(obj, 100);
  obj.icntl[block_size] = 4;
  obj.icntl[part_type] = 3;
  obj.dcntl[part_imbalance] = 2.5;

  EXPECT_NO_THROW(obj(-1));
  EXPECT_NO_THROW(obj(6));
}

TEST_F (AbcdTest, BlockCG_PaToH_5) 
{
  init_2d_lap(obj, 100);
  obj.icntl[block_size] = 4;
  obj.icntl[part_type] = 3;
  // this should make the solver throw an exception as it might produce an empty partition
  obj.dcntl[part_imbalance] = 5;

  EXPECT_NO_THROW(obj(-1));
  EXPECT_ANY_THROW(obj(6));
}

TEST_F (AbcdTest, CijAugment)
{
  init_2d_lap(obj, 100);
  obj.icntl[aug_type] = 1;

  EXPECT_NO_THROW(obj(-1));
  EXPECT_NO_THROW(obj(6));
}

TEST_F (AbcdTest, AijAugment)
{
  init_2d_lap(obj, 100);
  obj.icntl[aug_type] = 2;
  obj.icntl[part_type] = 3;
  obj.icntl[verbose_level] = 2;

  EXPECT_NO_THROW(obj(-1));
  EXPECT_NO_THROW(obj(6));
}

TEST_F (AbcdTest, AijAugment_SmallChuncks)
{
  init_2d_lap(obj, 10);
  obj.icntl[aug_type] = 2;
  obj.icntl[part_type] = 3;
  obj.icntl[aug_blocking] = 1;

  EXPECT_NO_THROW(obj(-1));
  EXPECT_NO_THROW(obj(6));
}

int main(int argc, char **argv) {
  // Equivalent to MPI_Initialize
  mpi::environment env(argc, argv);

  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

void init_2d_lap(abcd &obj, int mesh_size)
{
  mpi::communicator world;

  if(obj.comm.rank() != 0) return;

  obj.m = mesh_size*mesh_size; // number of rows
  obj.n = obj.m; // number of columns
  obj.nz = 3*obj.m - 2*mesh_size; // number of nnz in the lower-triangular part
  obj.sym = true;

  // allocate the arrays
  obj.irn = new int[obj.nz];
  obj.jcn = new int[obj.nz];
  obj.val = new double[obj.nz];

  init_2d_lap(obj.m, obj.n, obj.nz, obj.irn, obj.jcn, obj.val, mesh_size);

  // set the rhs
  obj.rhs = new double[obj.m];
  for (size_t i = 0; i < obj.m; i++) {
    obj.rhs[i] = ((double) i + 1)/obj.m;
  }
}

void init_2d_lap(int m, int n, int nz, int *irn, int *jcn, double *val, int mesh_size)
{
  int pos;
  int i;

  pos = 0;
  for (i = 1; i <= m; i++) {

    // the diagonal
    irn[pos] = i;
    jcn[pos] = i;
    val[pos] = -4.0;

    pos++;

    if (i % mesh_size != 0) {
      // the lower-triangular part
      irn[pos] = i + 1;
      jcn[pos] = i ;
      val[pos] = 1.0;
      pos++;
    } 

    if (i + mesh_size > m) continue;
    irn[pos] = i + mesh_size ;
    jcn[pos] = i ;
    val[pos] = 1.0;
    pos++;
  }
}
