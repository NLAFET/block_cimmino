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
 * \file solve.cpp
 * \brief Implementation of the solving using the Augmented Block Cimmino method
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include <abcd.h>
#include <iostream>
#include <fstream>

/*!
 *  \brief Build and solve S to compute the solution X=Abar^+b+W^+f
 *
 *  Build and solve S to compute the solution X=Abar^+b+W^+f
 *  with Abar^+=sum Abari^+bi
 *       W^+=(I-P)Y^TS^{-1}Ysum Abari^+bi
 *
 *  \param b: the right hand side
 */
void abcd::solveABCD ( MV_ColMat_double &b )
{
    /* Compute Infinite norm of the RHS */
    {
        // get first column of the local RHS
        MV_ColMat_double u(m, nrhs, 0);
        u = b(MV_VecIndex(0, b.dim(0)-1), MV_VecIndex(0,nrhs-1));

        // Parallel computation of the Inf-norm of the RHS
        nrmB = std::vector<double>(nrhs, 0);
        for(int j = 0; j < nrhs; ++j) {
            VECTOR_double u_j = u(j);
            double lnrmBs = infNorm(u_j);

            // Sync B norm :
            mpi::all_reduce(inter_comm, &lnrmBs, 1,  &nrmB[0] + j, mpi::maximum<double>());
        }
    }

    /* Compute w=A^+b */
    double t, tto = MPI_Wtime();
    MV_ColMat_double w;
    if(inter_comm.rank() == 0){
        LINFO << "*----------------------------------*";
        LINFO << "> Computing w = A^+b                ";
    }

    t = MPI_Wtime();
#ifdef WIP
    // If the augmentation is not filtered, sum of projections
    if(dcntl[Controls::aug_filter] == 0){
#endif //WIP

        w = sumProject(1e0, b, 0e0, Xk);

#ifdef WIP
    // else if filtered ABCD, use iterative BCG
    } else {
        bcg(b);
        w = Xk;
    }
#endif //WIP

    if(inter_comm.rank() == 0){
        LINFO << "> Time to compute w = A^+ b: " << setprecision(2) << MPI_Wtime() - t;
    }

    /* Compute the RHS second part f=-Yw */
    if(inter_comm.rank() == 0){
        LINFO << "*----------------------------------*";
        LINFO << "> Setting f = - Y w                 ";
    }

    // local f=-Abar^+b restricted to the augmented variables (last size_c dims)
    t = MPI_Wtime();
    //TODO: change 1 to nrhs
    MV_ColMat_double f(size_c, 1, 0); // local RHS
    for(std::map<int,int>::iterator it = glob_to_local.begin();
        it != glob_to_local.end(); ++it){
        if(it->first >= n_o){
            f(it->first - n_o, 0) = -1 * w(it->second, 0);
        }
    }

    // Reduce all local parts to get the global f of size size_c
    {
        double *f_ptr = f.ptr(); // pointer to local RHS
        MV_ColMat_double ff(size_c, 1, 0); // global RHS
        double *f_o = ff.ptr(); // pointer to global RHS
        mpi::all_reduce(inter_comm, f_ptr, size_c, f_o, or_bin);
        f = ff; // merge local to global
    }
    if(inter_comm.rank() == 0){
        LINFO << "> Time to centralize f : " << setprecision(2) << MPI_Wtime() - t;
    }

    /* Solve Sz=f */
    if(inter_comm.rank() == 0){
        LINFO << "*----------------------------------*";
        LINFO << "> Solving Sz = f                    ";
        LINFO << "*                                  *";
    }

    t = MPI_Wtime();
#ifdef WIP
    // if iterative solution of S, launch PCG(S)
    if(icntl[Controls::aug_iterative] != 0){
        if(inter_comm.rank() == 0)
            LINFO << "* ITERATIVELY                      *";

        VECTOR_double f0 = f.data();
        f0 = pcgS(f0);
        f.setData(f0);

    // else direct solve of S with MUMPS
    } else {

#endif //WIP

        f = solveS(f);

#ifdef WIP
    }
#endif //WIP

    if(IRANK == 0)
        LINFO << "| Time to solve Sz = f: " << setprecision(2) << MPI_Wtime() - t;

    /* Compute zz = (I-P)Y^Tz */
    t = MPI_Wtime();

    if(inter_comm.rank() == 0){
        LINFO << "*----------------------------------*";
        LINFO << "                                T   ";
        LINFO << "> [->] Computing zz = (I - P) Y  z  ";
    }

    // Xk=f
    Xk = MV_ColMat_double(n, 1, 0);
    MV_ColMat_double zrhs(m, 1, 0);

    for( int i = 0; i < size_c; i++){

        std::map<int,int>::iterator iti = glob_to_local.find(n_o + i);

        if(iti!=glob_to_local.end()){
            Xk(glob_to_local[n_o + i], 0) = f(i, 0);
        } else {
            continue;
        }
    }

#ifdef WIP
    if(dcntl[Controls::aug_filter] == 0){
#endif //WIP
        // f = f - sum Ai^+Aif
        f = Xk - sumProject(0e0, b, 1e0, Xk);

#ifdef WIP
    // Filters the augmentation if needed
    } else {
        use_xk = false;
        f = MV_ColMat_double(n, 1, 0);

        // ????
        if(!use_xk){
            int st = 0;
            for(int p = 0; p < nb_local_parts; p++){
                MV_ColMat_double sp =  spsmv(partitions[p], local_column_index[p], Xk);
                int pos = 0;
                for(int k = st; k < st + partitions[p].dim(0); k++){
                    zrhs(k, 0) = sp(pos, 0);
                    pos++;
                }
                st += partitions[p].dim(0);
            }

            for(int i = 0; i < f.dim(0); i++)
                f(i, 0) = Xk(i, 0);
        }

        bcg(zrhs);

        if(use_xk)
            f = Xk;
        else
            f = f - Xk;
    }
#endif // WIP
    if(IRANK == 0)
        LINFO << "> Took: " << setprecision(2) << MPI_Wtime() - t;

    /*
     * The final solution (distributed)
     * w = \Abar^+ b
     * f = \Wbar^+ b
     */
    Xk = w + f;

    if(IRANK == 0){
        LINFO << "*----------------------------------*";
        LINFO << "> Total time to solve: " << setprecision(2) << MPI_Wtime() - tto;
        LINFO << "*----------------------------------*";
        LINFO << "";
    }

    // Compute backward error
    compute_rho(Xk, b);

    // Centralize solution
    if(IRANK == 0) {
        solution = MV_ColMat_double(n_o, nrhs, 0);
        sol = solution.ptr();
    }
    centralizeVector(sol, n_o, nrhs, Xk.ptr(), n, nrhs, glob_to_local_ind, &dcol_[0]);

}		/* -----  end of function abcd::solveABCD  ----- */
