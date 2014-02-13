#include "vect_utils.h"
#include <climits>
#include <iostream>
#include <assert.h>
using namespace std;

//! \brief Returns the index of non-null columns
//! @param[in] col_ptr Compressed column array
//! @param[in] _size   Size - 1 of the @col_ptr array
//! \return Indices of non-null columns
std::vector<int> getColumnIndex(int *col_ptr, int _size)
{
    std::vector<int> column_index;
    int *index1 = col_ptr,
        *index2 = col_ptr + 1,
        *index3 = col_ptr + _size + 1,
        j = 0;
        
    while(index2 != index3)
    {
        if(*index1 != *index2)
        {
            column_index.push_back(j);
            index1 = index2;
        }
        index2++;
        j++;
    }
    return column_index;
}

//! \brief Merge sorted vectors and keep the sorting
//! @param[in] vectors the vectors to be merged
//! \return The merge
std::vector<int> mergeSortedVectors(std::vector<std::vector<int> > &vectors)
{
    if(vectors.size() == 1) return vectors[0];
    std::vector<int> merge;
    
    std::map<int,std::vector<int>::iterator> its;
    std::map<int,std::vector<int>::iterator> eds;
    for(int i = 0; i < vectors.size(); i++)
    {
        its[i] = vectors[i].begin();
        eds[i] = vectors[i].end();
    }

    int min;
    std::vector<int> to_delete;
    while(its.size() != 0)
    {
        min = INT_MAX;
        std::vector<int>::iterator min_it;
        int idx = -1;

        for(std::map<int,std::vector<int>::iterator>::iterator it = its.begin();
                it != its.end(); it++)
        {
            if(*it->second <= min){
                min_it = it->second;
                min = *it->second;
                idx = it->first;
            }
        }
        for(std::map<int,std::vector<int>::iterator>::iterator it = its.begin();
                it != its.end(); it++)
        {
                if(*it->second == min && it->second + 1 == eds[it->first])
                {
                    to_delete.push_back(it->first);
                }
        }

#ifdef DEBUG_ASSERTS
        assert(eds[idx] - its[idx] >= 0);
#endif

        for(std::map<int,std::vector<int>::iterator>::iterator it = its.begin();
                it != its.end(); it++)
        {
            if(*it->second == min && it->second != min_it){
                it->second++;
            } else if(it->second == eds[it->first]){
                to_delete.push_back(it->first);
            }
        }

        merge.push_back(*min_it);
        //++min_it;
        its[idx]++;
        if(to_delete.size() != 0)
            for(std::vector<int>::iterator it = to_delete.begin();
                    it!= to_delete.end(); it++)
                its.erase(*it);
        to_delete.clear();
    }
    return merge;
}

//! \brief Returns the corresponding indices to an intersection
//! @param[in] v1 first vector 
//! @param[in] v2 second vector 
//! \return The indices of intersection
std::pair<std::vector<int>, std::vector<int> >
    getIntersectionIndices(std::vector<int> &v1, std::vector<int> &v2)
{

    std::vector<int> inter1;
    std::vector<int> inter2;
    std::vector<int>::iterator it1 = v1.begin();
    std::vector<int>::iterator it2 = v2.begin();

    inter1.reserve(v1.size());
    inter2.reserve(v2.size());

    int id1 = 0, id2 = 0;

    while(it1 != v1.end() && it2 != v2.end()) {
        if(*it1 < *it2) {
            ++it1;
            ++id1;
        } else {
            if(!(*it2 < *it1)) {
                inter1.push_back( id1);
                inter2.push_back( id2);

                ++it1; ++id1;
            }
            ++id2;
            ++it2;
        }
    }
    std::pair<std::vector<int>, std::vector<int> > intersection(inter1, inter2);

    return intersection;
}
