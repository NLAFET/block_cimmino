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

    int job = 0;
    do{
        mpi::broadcast(intra_comm, job, 0);
        
        if(job == -1) break;
        if(job == 1){
            mumps.job = 3;
            dmumps_c(&mumps);
        } else if (job == 2) {
            mu.sym = 2;
            mu.par = 1;
            mu.job = -1;

            mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) world);
            dmumps_c(&mu);

            mu.icntl[0] = -1;
            mu.icntl[1] = -1;
            mu.icntl[2] = -1;

            mu.job = 1;
            dmumps_c(&mu);
            mu.job = 2;
            dmumps_c(&mu);
            mu.job = 3;
            dmumps_c(&mu);
        }
    }while(true);

}		/* -----  end of function abcd::waitForSolve()  ----- */



