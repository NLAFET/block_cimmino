#include <abcd.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::waitForSolve
 *  Description:  
 * =====================================================================================
 */
    void
abcd::waitForSolve()
{
    DMUMPS_STRUC_C mu;
    mpi::communicator world;

    vector<int> vrows, vcols;
    vector<double> vvals;

    int job = 0;
    do{
        mpi::broadcast(intra_comm, job, 0);
        
        //we are finished here!
        if(job == -1) break;
        // do the solve, in the (simple or not) sumproject
        if(job == 1){
            mumps.job = 3;
            dmumps_c(&mumps);
        // handle the matrix S
        } else if (job == 2) {

            mu.sym = 2;
            mu.par = 1;
            mu.job = -1;

            mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) world);
            dmumps_c(&mu);

            mu.icntl[0] = -1;
            mu.icntl[1] = -1;
            mu.icntl[2] = -1;

            mu.nz_loc = vrows.size();
            mu.irn_loc = &vrows[0];
            mu.jcn_loc = &vcols[0];
            mu.a_loc = &vvals[0];

            mu.job = 1;
            dmumps_c(&mu);
            mu.job = 2;
            dmumps_c(&mu);
            mu.job = 3;
            dmumps_c(&mu);
        // do the solve with a distributed solution output, when building S
        } else if (job == 3) {
            int s;
            int start_c;
            // get some necessary arrays and data
            //
            mpi::broadcast(intra_comm, s, 0);
            std::vector<int> loc_cols, mycols;
            mpi::broadcast(intra_comm, start_c, 0);

            if(start_c == -1)
                continue;

            mpi::broadcast(intra_comm, column_index, 0);
            mpi::broadcast(intra_comm, loc_cols, 0);
            mpi::broadcast(intra_comm, n_o, 0);
            mpi::broadcast(intra_comm, mycols, 0);

            mumps.lsol_loc = mumps.getInfo(23);
            mumps.sol_loc = new double[mumps.lsol_loc * s];
            mumps.isol_loc = new int[mumps.lsol_loc];

            mumps.job = 3;
            dmumps_c(&mumps);

            int part = 0;
            int x_pos = 0;

            x_pos = start_c;
            int end_c = column_index[part].size();
            int i_loc = 0;
            vector<int> target;
            vector<int> target_idx;
            target.reserve(size_c);
            target_idx.reserve(size_c);

            while (i_loc < mumps.lsol_loc){
                int isol = mumps.isol_loc[i_loc] - 1;
                if (isol >= start_c && isol < end_c){
                    int ci = column_index[part][isol] - n_o; 
                    target.push_back(i_loc);
                    target_idx.push_back(ci);
                }
                i_loc++;
            }

            for (int j = 0; j < s; j++) {
                int col = mycols[j];
                for (int i = 0; i < (int)target_idx.size(); i++) {
                    if(target_idx[i] < col)
                        continue;
                    else if ( target_idx[i] == col)
                        vvals.push_back(0.5 - mumps.sol_loc[target[i] + j * mumps.lsol_loc]);
                    else
                        vvals.push_back(-mumps.sol_loc[target[i] + j * mumps.lsol_loc]);
                    vrows.push_back(target_idx[i] + 1);
                    vcols.push_back(col + 1);

                }
            }

            delete[] mumps.isol_loc;
            delete[] mumps.sol_loc;
        }
    }while(true);

}		/* -----  end of function abcd::waitForSolve()  ----- */



