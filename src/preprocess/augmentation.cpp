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
 * \file augmentation.cpp
 * \brief Implementation of the Gateway function to call the actual augmentation
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include "abcd.h"

/*!
 *  \brief Augments the matrix and build the C part in [A C]
 *
 *  Augment the matrix M by building the augmentation part C. This method
 *  calls the chosen specific implementation:
 *    - aug_type=0 => no augmentation
 *    - aug_type=1 => Aij augmentation
 *    - aug_type=2 => Cij augmentation
 *
 *  \param M: A vector of the compressed partitions (splib format) to augment
 *
 */
void abcd::augmentMatrix (std::vector<CompCol_Mat_double> &M)
{
    stC = std::vector<int>(M.size(), -1);

std::cout << "AUGMENTATION TYPE: " << icntl[Controls::aug_type] << "\n";
    switch(icntl[Controls::aug_type]) {
        // No augmentation
        case 0: {
            return;
        // C_ij/-I augmentation
        } case 1: {
            cijAugmentMatrix(M);
            break;
        // A_ij/-A_ji augmentation
        } case 2: {
            aijAugmentMatrix(M);
            break;
        } default: {
            info[Controls::status] = -8;
            mpi::broadcast(comm, info[Controls::status], 0);
            throw std::runtime_error("Unkown augmentation scheme.");
        }
    }
}               /* -----  end of function abcd::augmentMatrix  ----- */
