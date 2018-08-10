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

#include <abcd.h>
#include <mumps.h>

/*!
 *  \brief create the augmented systems for the local partitions
 *
 *  Create block diagonal matrix with the augmented systems of all local partitions
 *
 *  \param n_aug: size of the block diagonal augmented system
 *  \param nz_aug: number of non-zeros of the block diagonal augmented system
 *  \param irn_aug: array of rows of the block diagonal augmented system
 *  \param jcn_aug: array of columns of the block diagonal augmented system
 *  \param val_aug: array of values of the block diagonal augmented system
 *
 */
void abcd::createAugmentedSystems(int &n_aug, int &nz_aug,
        std::vector<int> &irn_aug, std::vector<int> &jcn_aug, std::vector<double> &val_aug)
{
    // Size of augmented system after gathering partitions
    m_n = 0; // number of columns in the augmented system
    m_nz = 0; // number of non-zeros in the augmented system
    for(int j = 0; j < nb_local_parts; j++) {
        m_n += partitions[j].dim(0) + partitions[j].dim(1);
        m_nz += partitions[j].dim(1) + partitions[j].NumNonzeros();
    }

    // Allocate the data for the augmented system
    n_aug = m_n;
    nz_aug = m_nz;
    irn_aug.resize(m_nz);
    jcn_aug.resize(m_nz);
    val_aug.resize(m_nz);

    // Use Fortran array (MUMPS) => start from 1
    int i_pos = 1;
    int j_pos = 1;
    int st = 0;

    // Build the augmented system
    for(int p = 0; p < nb_local_parts; ++p) {
        // fill the identity
        for(int i = 0; i < partitions[p].dim(1); ++i) {
            irn_aug[st + i] = i_pos + i;
            jcn_aug[st + i] = j_pos + i;
            val_aug[st + i] = dcntl[Controls::alpha];
        }

        // we get down by nb_cols
        i_pos += partitions[p].dim(1);
        // we added nb_cols elements
        st += partitions[p].dim(1);

        // Add partition in lower triangular part (symmetric augmented systems)
        for(int k = 0; k < partitions[p].dim(0); ++k) {
            for(int j = partitions[p].row_ptr(k); j < partitions[p].row_ptr(k + 1); ++j) {
                irn_aug[st] = i_pos + k;
                jcn_aug[st] = j_pos + partitions[p].col_ind(j);
                val_aug[st] = partitions[p].val(j);

                st++;
            }
        }

        // shift to build augmented system of next partition
        i_pos += partitions[p].dim(0);
        j_pos += partitions[p].dim(1) + partitions[p].dim(0);

    }
}               /* -----  end of function abcd::createAugmentedSystems  ----- */
