#include <abcd.h>
#include <mumps.h>

void abcd::createAugmentedSystems(int &n_aug, int &nz_aug,
        std::vector<int> &irn_aug, std::vector<int> &jcn_aug, std::vector<double> &val_aug)
{
    m_n = 0;
    m_nz = 0;

    for(int j = 0; j < nb_local_parts; j++) {
        m_n += partitions[j].dim(0) + partitions[j].dim(1);
        m_nz += partitions[j].dim(1) + partitions[j].NumNonzeros();
    }

    // Allocate the data for mu
    n_aug = m_n;
    nz_aug = m_nz;
    irn_aug.resize(m_nz);
    jcn_aug.resize(m_nz);
    val_aug.resize(m_nz);

    // Use Fortran array => start from 1
    int i_pos = 1;
    int j_pos = 1;
    int st = 0;

    for(int p = 0; p < nb_local_parts; ++p) {

        // fill the identity
        for(int i = 0; i < partitions[p].dim(1); ++i) {
            irn_aug[st + i] = i_pos + i;
            jcn_aug[st + i] = j_pos + i;
            val_aug[st + i] = 1;
        }

        // we get down by nb_cols
        i_pos += partitions[p].dim(1);
        // we added nb_cols elements
        st += partitions[p].dim(1);

        for(int k = 0; k < partitions[p].dim(0); ++k) {
            for(int j = partitions[p].row_ptr(k); j < partitions[p].row_ptr(k + 1); ++j) {
                irn_aug[st] = i_pos + k;
                jcn_aug[st] = j_pos + partitions[p].col_ind(j);
                val_aug[st] = partitions[p].val(j);

                st++;
            }
        }

        i_pos += partitions[p].dim(0);
        j_pos += partitions[p].dim(1) + partitions[p].dim(0);

    }
}
