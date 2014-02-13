#include <abcd.h>
#include <mumps.h>

void abcd::createAugmentedSystems(MUMPS &mu)
{
    m_n = 0;
    m_nz = 0;

    for(int j = 0; j < nb_local_parts; j++) {
        m_n += partitions[j].dim(0) + partitions[j].dim(1);
        m_nz += partitions[j].dim(1) + partitions[j].NumNonzeros();
    }

    // Allocate the data for mu
    mu.n = m_n;
    mu.nz = m_nz;
    mu.irn = new int[m_nz];
    mu.jcn = new int[m_nz];
    mu.a = new double[m_nz];

    // Use Fortran array => start from 1
    int i_pos = 1;
    int j_pos = 1;
    int st = 0;

    for(int p = 0; p < nb_local_parts; p++) {

        // fill the identity
        for(int i = 0; i < partitions[p].dim(1); i++) {
            mu.irn[st + i] = i_pos + i;
            mu.jcn[st + i] = j_pos + i;
            mu.a[st + i] = 1;
        }

        // we get down by nb_cols
        i_pos += partitions[p].dim(1);
        // we added nb_cols elements
        st += partitions[p].dim(1);

        for(int k = 0; k < partitions[p].dim(0); k++) {
            for(int j = partitions[p].row_ptr(k); j < partitions[p].row_ptr(k + 1); j++) {
                mu.irn[st] = i_pos + k;
                mu.jcn[st] = j_pos + partitions[p].col_ind(j);
                mu.a[st] = partitions[p].val(j);

                st++;
            }
        }

        i_pos += partitions[p].dim(0);
        j_pos += partitions[p].dim(1) + partitions[p].dim(0);

    }
}
