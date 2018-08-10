// Copyright Institut National Polytechnique de Toulouse (2014) 
// Contributor(s) :
// M. Zenadi <mzenadi@enseeiht.fr>
// D. Ruiz <ruiz@enseeiht.fr>
// R. Guivarch <guivarch@enseeiht.fr>

// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html"

// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 

// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 

// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.

/*!
 * \file compress_columns.cpp
 * \brief Implementation of some functions to compress the columns of a matrix
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

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
}               /* -----  end of function getColumnIndex  ----- */

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
    while(!its.empty())
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
}               /* -----  end of function mergeSortedVectors  ----- */

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
}               /* -----  end of function getIntersectionIndices  ----- */
