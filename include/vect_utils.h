#ifndef VECT_UTILS_HXX_
#define VECT_UTILS_HXX_

#include<vector>
#include<map>

std::vector<int> getColumnIndex(int *col_ptr, int _size);

std::vector<int> mergeSortedVectors(std::vector<std::vector<int> > &vectors);

std::pair<std::vector<int>, std::vector<int> >
    getIntersectionIndices(std::vector<int> &v1, std::vector<int> &v2);


#endif // MAT_UTILS_HXX_
