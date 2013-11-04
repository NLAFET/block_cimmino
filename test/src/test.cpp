#include "abcd.h"
#include "mat_utils.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
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
