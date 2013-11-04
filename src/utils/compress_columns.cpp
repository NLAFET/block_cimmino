#include "vect_utils.h"
#include <climits>
#include <iostream>

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
    std::vector<int> merge;
    if(vectors.size() == 1) return vectors[0];
    
    std::map<int,std::vector<int>::iterator> its;
    std::map<int,std::vector<int>::iterator> eds;
    for(int i = 0; i < vectors.size(); i++)
    {
        its[i] = vectors[i].begin();
        eds[i] = vectors[i].end();
    }

    int min, min_index, to_delete = -1;
    while(its.size() != 0)
    {
        min = INT_MAX;
        to_delete = -1;
        min_index = -1;
        for(std::map<int,std::vector<int>::iterator>::iterator it = its.begin();
                it != its.end(); it++)
        {
            if(*it->second < min){
                min_index = it->first;
                min = *it->second;

                if(it->second + 1 == eds[it->first])
                {
                    to_delete = it->first;
                }
            } else if(*it->second == min) it->second++;
        }

        merge.push_back(*its[min_index]);
        its[min_index]++;
        if(to_delete != -1) its.erase(to_delete);
    }
    return merge;
}
