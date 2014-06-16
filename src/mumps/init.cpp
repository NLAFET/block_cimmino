#include <mumps.h>
#include <abcd.h>

void abcd::initializeMumps(MUMPS &mu)
{
    initializeMumps(mu, false);
}

void abcd::initializeMumps(MUMPS &mu, bool local)
{
    // The first run of MUMPS is local to CG-masters
    std::vector<int> r;
    if(local) {
        r.push_back(inter_comm.rank());
        mpi::group grp = inter_comm.group().include(r.begin(), r.end());
        intra_comm = mpi::communicator(inter_comm, grp);
    } else {
        if(instance_type == 0) {
            r.push_back(comm.rank());
        } else {
            r.push_back(my_master);
        }
        std::copy(my_slaves.begin(), my_slaves.end(), std::back_inserter(r));

        mpi::group grp = comm.group().include(r.begin(), r.end());
        intra_comm = mpi::communicator(comm, grp);
    }

    mu.sym = 2;
    mu.par = 1;
    mu.job = -1;
    mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) intra_comm);

    dmumps_c(&mu);
    if(mu.getInfo(1) != 0) throw mu.getInfo(1);

    mu.initialized = true;

    mu.setIcntl(1, -1);
    mu.setIcntl(2, -1);
    mu.setIcntl(3, -1);

    mu.setIcntl(6, 5);
    mu.setIcntl(7, 5);
    mu.setIcntl(8, -2);
    mu.setIcntl(12, 2);
    mu.setIcntl(14, 90);
    if(inter_comm.rank() == 0) strcpy(mu.write_problem,"/tmp/test.mm");
}
