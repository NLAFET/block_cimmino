#include "abcd.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "vect_utils.h"

using ::testing::AtLeast;
using ::testing::Return;
using ::testing::Eq;

TEST (getColumnIndex, MatrixWithSomeEmptyCols) { 
    int sz = 10;
    int col_ptr[11] = {0, 2, 2, 3, 3, 4, 5, 5, 7, 8, 9};
    std::vector<int> ci;
    ci.reserve(sz);
    
    int j = 0;
    for(int i = 1; i <= sz; i++) {
        if(col_ptr[i] != col_ptr[i - 1]) ci.push_back(j);
        j++;
    }
    // ci should be equal to {0, 2, 4, 5, 7, 8, 9}!
    EXPECT_THAT(getColumnIndex(col_ptr, sz),Eq(ci));
}

TEST (getColumnIndex, MatrixWithNoEmptyCols) { 
    int sz = 10;
    int col_ptr[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<int> ci;
    ci.reserve(sz);
    
    int j = 0;
    for(int i = 1; i <= sz; i++) {
        if(col_ptr[i] != col_ptr[i - 1]) ci.push_back(j);
        j++;
    }
    // ci should be equal to {0,1,2,3,4,5,6,7,8,9}!
    EXPECT_THAT(getColumnIndex(col_ptr, sz),Eq(ci));
}

TEST (getColumnIndex, MatrixWithOnlyEmptyCols) { 
    int sz = 10;
    int col_ptr[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> ci;
    ci.reserve(sz);
    
    int j = 0;
    for(int i = 1; i <= sz; i++) {
        if(col_ptr[i] != col_ptr[i - 1]) ci.push_back(j);
        j++;
    }
    // ci should be equal to {}!
    EXPECT_THAT(getColumnIndex(col_ptr, sz),Eq(ci));
}

TEST (mergeSortedVectors, SingleVector) { 
    int v1[3] = {1, 3, 5};
    std::vector<int> vv1(&v1[0], &v1[0]+3);
    std::vector<std::vector<int> > v;
    v.push_back(vv1);

    EXPECT_THAT(mergeSortedVectors(v), Eq(vv1));
}

TEST (mergeSortedVectors, TestMerge) { 
    int v1[3] = {1, 2, 5};
    std::vector<int> vv1(&v1[0], &v1[0]+3);
    int v2[4] = {2, 4, 6, 8};
    std::vector<int> vv2(&v2[0], &v2[0]+4);
    int v3[7] = {1, 2, 3, 4, 6, 8, 9};
    std::vector<int> vv3(&v3[0], &v3[0]+7);
    
    std::vector<std::vector<int> > v;
    v.push_back(vv1); v.push_back(vv2); v.push_back(vv3);
    
    int r[8] = {1, 2, 3, 4, 5, 6, 8, 9};
    std::vector<int> re(&r[0], &r[0]+8);
    
    EXPECT_THAT(mergeSortedVectors(v), Eq(re));
}

TEST (intersectAllToAll, TestIntersect) { 
    int v1[4] = {0, 2, 3, 5};
    std::vector<int> vv1(v1, v1+4);
    int v2[5] = {0, 1, 2, 4, 8};
    std::vector<int> vv2(v2, v2+5);
    int v3[4] = {4, 6, 8, 9};
    std::vector<int> vv3(v3, v3+4);

    int v12[2] = {0, 1};
    std::vector<int> vv12(v12, v12+2);
    int v21[2] = {0, 2};
    std::vector<int> vv21(v21, v21+2);

    std::vector<int> vv13;
    std::vector<int> vv31;

    int v23[2] = {3, 4};
    std::vector<int> vv23(v23, v23+2);
    int v32[2] = {0, 2};
    std::vector<int> vv32(v32, v32+2);

    std::pair<std::vector<int>, std::vector<int> > p12(vv12, vv21), p13(vv13, vv31), p23(vv23, vv32);

    EXPECT_THAT(getIntersectionIndices(vv1, vv2), Eq(p12));
    EXPECT_THAT(getIntersectionIndices(vv1, vv3), Eq(p13));
    EXPECT_THAT(getIntersectionIndices(vv2, vv3), Eq(p23));
}
