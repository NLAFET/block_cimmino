#include <sparse/sparse.hpp>
#include <sparse/coo.hpp>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <iostream>

using ::testing::AtLeast;
using ::testing::Return;
using ::testing::Eq;

TEST (CooMatrix_Constructor, CoordinateConstructor)
{
    int m = 5;
    int n = 5;
    int nz = 8;
    int irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    int jcn[8] = {0, 1, 2, 4, 3, 4, 2, 1};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CooMatrix<> cm(m, n, nz, irn, jcn, val);
}

TEST (CooMatrix_Constructor, CoordinateConstructorCustomType)
{
    size_t m = 5;
    size_t n = 5;
    size_t nz = 8;
    size_t irn[8] = {0, 1, 2, 2, 3, 3, 4, 4};
    size_t jcn[8] = {0, 1, 2, 4, 3, 4, 2, 1};
    double val[8] = {1, 1, 1, 3, 1, 2, 3, 2};

    CooMatrix<double,size_t,0> cm(m, n, nz, irn, jcn, val);
}
