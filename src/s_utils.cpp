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
 * \file s_utils.cpp
 * \brief Implementation of utils to build/solve S directly or iteratively, preconditioned with M or not
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include <abcd.h>
#include <iostream>
#include <fstream>

#ifdef WIP
/*!
 *  \brief Compute the product S*V
 *
 *  Compute the product SV=(I-sum Abari^+)V
 *
 *  \param V: Vector which product with S must be computed
 *
 */
    MV_ColMat_double
abcd::prodSv ( MV_ColMat_double &V )
{
    MV_ColMat_double W(n, V.dim(1), 0); // local V
    MV_ColMat_double b; // just for decoration!
    MV_ColMat_double R(V.dim(0), V.dim(1), 0); // SV

    // W=V in local indexing
    for(std::vector<int>::iterator iti = st_c_part_it; iti != glob_to_local_ind.end(); ++iti){
            for(int j = 0; j < V.dim(1); j++)
                W(iti - glob_to_local_ind.begin() , j) = V(*iti - n_o, j);
    }

    // W=(I-sum Abari^+)W
    W = W - sumProject(0e0, b, 1e0, W);

    // R=W in global indexing
    for(std::vector<int>::iterator iti = st_c_part_it; iti != glob_to_local_ind.end(); ++iti){
            for(int j = 0; j < V.dim(1); j++)
                R(*iti - n_o, j) = W(iti - glob_to_local_ind.begin() , j);

    }

    return R;
}       /* -----  end of function abcd::prodSv  ----- */

/*!
 *  \brief Build the preconditioning matrix M and its corresponding MUMPS
 *
 *  Build the preconditioning matrix M (restriction on rows and columns) and 
 *  its corresponding MUMPS
 *
 */
MUMPS abcd::buildM (  )
{
    if(inter_comm.rank() == 0){
        clog << "* Building the preconditioner      *" << endl;
        clog << " Size is : " << selected_S_columns.size() << endl;

        //ofstream f;
        //f.open("/tmp/selected");
        //for(int i = 0; i < selected_S_columns.size(); i++){
            //f << selected_S_columns[i] << endl;
        //}
        //f.close();
    }


    double t = MPI_Wtime();
    //t = MPI_Wtime();
    //S = buildS(skipped_S_columns);
    //clog << " Time to build S : " << MPI_Wtime() - t << endl;

    // Build S on selected columns
    Coord_Mat_double M = buildS(selected_S_columns);
    //clog << " Time to build M : " << MPI_Wtime() - t << endl;

    //if(IRANK == 0) clog << "* -----------  DONE  ------------- * " << endl;

    std::vector<int> dropped;

    {
        // arrays of rows, columns and values for M
        std::vector<int> mr, mc;
        std::vector<double> mv;

        std::vector<int>::iterator it;

        // M is 1 on skipped columns of S
        for(size_t i = 0; i < skipped_S_columns.size(); i++){
            std::map<int,int>::iterator iti = glob_to_local.find(n_o + skipped_S_columns[i]);
            if(iti == glob_to_local.end()) continue;

            int ro = skipped_S_columns[i];

            mr.push_back(ro);
            mc.push_back(ro);
            //mv.push_back(S(ro,ro));
            mv.push_back(1);
        }

        // ?????
        for(int i = 0; i < M.NumNonzeros(); i++){
            if(std::find(skipped_S_columns.begin(), skipped_S_columns.end(), M.row_ind(i))
                    !=skipped_S_columns.end()) continue;
            if(M.val(i) == 0) continue;
            mr.push_back(M.row_ind(i));
            mc.push_back(M.col_ind(i));
            mv.push_back(M.val(i));
        }

        M = Coord_Mat_double(size_c, size_c, mv.size(), &mv[0], &mr[0], &mc[0]);
    }

    /*-----------------------------------------------------------------------------
     *  MUMPS initialization
     *-----------------------------------------------------------------------------*/
    MUMPS mu;
    mu.sym = 1;
    mu.par = 1;
    mu.job = -1;
    mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) inter_comm);

    dmumps_c(&mu);

    mu.icntl[0] = -1;
    mu.icntl[1] = -1;
    mu.icntl[2] = -1;

    if(inter_comm.rank() == 0){
        //strcpy(mu.write_problem, "/tmp/mmm.mtx");
        //mu.icntl[0] = 6;
        //mu.icntl[1] = 6;
        //mu.icntl[2] = 6;
        //mu.icntl[3] = 2;
    }

    mu.n = M.dim(0);

    // parallel analysis if the S is large enough
    //if(mu.n >= 200) {
        //mu.icntl[28 - 1] =  2;
    //}
    mu.icntl[8  - 1] =  7;
    mu.icntl[7  - 1] =  5;
    mu.icntl[14 - 1] =  90;

    if(inter_comm.size() == 1){
        mu.nz= M.NumNonzeros();
        mu.irn= M.rowind_ptr();
        mu.jcn= M.colind_ptr();
        mu.a= M.val_ptr();
        for(int i = 0; i < mu.nz; i++){
            mu.irn[i]++;
            mu.jcn[i]++;
        }
    } else {
        mu.icntl[18 - 1]= 3;
        mu.nz_loc = M.NumNonzeros();

        mu.irn_loc = M.rowind_ptr();
        mu.jcn_loc = M.colind_ptr();
        mu.a_loc = M.val_ptr();
        for(int i = 0; i < mu.nz_loc; i++){
            mu.irn_loc[i]++;
            mu.jcn_loc[i]++;
        }
    }

    /*-----------------------------------------------------------------------------
     *  MUMPS analysis
     *-----------------------------------------------------------------------------*/
    t = MPI_Wtime();

    IBARRIER;
    if(inter_comm.rank() == 0)
        clog << "* Starting Analysis of the precond *" << endl;

    mu(1)

    if(inter_comm.rank() == 0){
        clog << "*                                  *" << endl;
        clog << "  [--] T.Analyse M :   " << MPI_Wtime() - t << endl;
        clog << "*                                  *" << endl;
    }

    /*-----------------------------------------------------------------------------
     *  MUMPS factorization
     *-----------------------------------------------------------------------------*/
    t = MPI_Wtime();

    mu(2);

    if(inter_comm.rank() == 0){
        clog << "*                                  *" << endl;
        clog << "  [--] T.Factorize M : " << MPI_Wtime() - t << endl;
        clog << "  N, NZ              : " << mu.n << ", " << mu.nz << endl;
        clog << "*                                  *" << endl;
    }

    mu.icntl[0] = -1;
    mu.icntl[1] = -1;
    mu.icntl[2] = -1;

    if(mu.info[0] < 0) {
        clog << IRANK << " 1 : " << mu.info[0] <<  " 2 : " << mu.info[1] << endl;
        exit(0);
    }
    /*-----------------------------------------------------------------------------
     *  END MUMPS part
     *-----------------------------------------------------------------------------*/
    return mu;
}       /* -----  end of function abcd::buildM  ----- */


/*!
 *  \brief Solve MUMPS on the preconditioning matrix M: M^{-1}*ones()
 *
 *  Solve MUMPS on the preconditioning matrix M: M^{-1}*ones()
 *
 */
    VECTOR_double
abcd::solveM (MUMPS &mu, VECTOR_double &z )
{
    // Initialize MUMPS RHS with z (residual)
    if(IRANK == 0) {
        mu.rhs = new double[z.size()];
        for(int i = 0; i < z.size(); i++) mu.rhs[i] = z(i);
        mu.nrhs = 1;
    }

    // Launch MUMPS solve on !!!!!!!!!!!
    mu(3);

    // Broadcast solution M^{-1}*ones()
    if(IRANK == 0){
        VECTOR_double sol(mu.rhs, z.size());
        double *sol_ptr = sol.ptr();
        mpi::broadcast(inter_comm, sol_ptr, size_c, 0);
        delete[] mu.rhs;
        return sol;
    } else {
        VECTOR_double sol(z.size(), 0);
        double *sol_ptr = sol.ptr();
        mpi::broadcast(inter_comm, sol_ptr, size_c, 0);
        return sol;
    }

}       /* -----  end of function abcd::solveM  ----- */

/*!
 *  \brief Solve on S using PCG
 *
 *  Solve on S using PCG possibly preconditioned with a matrix M
 *
 */
    VECTOR_double
abcd::pcgS ( VECTOR_double &b )
{
    double resid; // scaled residual
//    double tol = 1e-5;
    double tol = dcntl[Controls::threshold]; // PCG threshold
    int max_iter = size_c; // maximum #iterations: size of the augmentation
    MUMPS mu;

    // If we want to solve S iteratively with preconditioning matrix M
    if(dcntl[Controls::aug_precond] != 0) mu = buildM();

    double t = MPI_Wtime();

    // PCG variables
    VECTOR_double p, z, q;
    VECTOR_double alpha(1), beta(1), rho(1), rho_1(1);

    // **************************************************
    // ITERATION k = 0                                 *
    // **************************************************
    // Compute Szk
    MV_ColMat_double zk(size_c, 1, 0);
    VECTOR_double x = zk.data();
    MV_ColMat_double Szk = prodSv(zk); // Sz_k

    // initialize residual vector and scaled_residual figure
    VECTOR_double r = b - Szk(0);
    double normb = norm(b); // norm of b
    if (normb == 0.0)
        normb = 1;
    resid = norm(r) / normb;
    // If scaled residual already under threshold, no need for PCG iterations
    if (resid <= tol) {
        tol = resid;
        max_iter = 0;
        return 0;
    }

    // **************************************************
    // ITERATIONs                                      *
    // **************************************************
    for (int i = 1; i <= max_iter; i++) {
        // z = M^{-1}R
        TIC;
        if(dcntl[Controls::aug_precond] != 0) {
            z = solveM(mu, r);
        } else {
            z = r;
        }
        //IFMASTER clog << " Time in solveM " << TOC << endl;

        TIC;
        // backward error r*z'
        rho(0) = dot(r, z);

        // Update direction p=z+beta*p (only z at the first iteration)
        if (i == 1)
            p = z;
        else {
            beta(0) = rho(0) / rho_1(0);
            p = z + beta(0) * p;
        }

        // Compute Szk
        MV_ColMat_double pv(p.ptr(), size_c, 1);
        //IFMASTER clog<< "Time in dot " << TOC << endl;
        TIC;
        Szk = prodSv(pv);
        //IFMASTER clog<< "Time in prodsv " << TOC << endl;
        TIC;

        // Reduce f over all masters
        double *f_ptr = Szk.ptr();
        MV_ColMat_double ff(size_c, 1, 0);
        double *f_o = ff.ptr();

        mpi::all_reduce(inter_comm, f_ptr, size_c, f_o, or_bin);

        // Q=SM^{-1}R
        q = ff(0);

        // alpha=RM^{-1}R / PQ
        alpha(0) = rho(0) / dot(p, q);

        //Update X and R
        x += alpha(0) * p;
        r -= alpha(0) * q;

        // COmpute scaled residual
        resid = norm(r) / normb;
        if (resid <= tol) {
            tol = resid;
            max_iter = i;
            if(IRANK == 0){
                clog << "Iteration to solve Sz = f " << i << " with a residual of " << resid << endl;
                clog << "Time in iterations : " << MPI_Wtime() - t << endl;
            }
            return x;
        }
        if(IRANK == 0)
            clog << "Iteration  " << i << " residual " << resid << "\r" << flush;

        //IFMASTER clog<< "Time in otherstuffs " << TOC << endl;

        // Save backward error
        rho_1(0) = rho(0);
    }

    // Display information on PCG for S
    tol = resid;
    if(IRANK == 0){
        clog << "Iteration to solve Sz = f " << max_iter << " with a residual of " << resid << endl;
        clog << "Time in iterations : " << MPI_Wtime() - t << endl;
    }

    return x;
}       /* -----  end of function abcd::pcgS  ----- */

#endif //WIP
