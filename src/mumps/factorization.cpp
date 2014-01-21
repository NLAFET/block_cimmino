#include<abcd.h>
#include<mumps.h>

void abcd::factorizeAugmentedSystems(MUMPS &mu)
{
    mpi::communicator world;
    mu.job = 2;

    double t = MPI_Wtime();
    if(inter_comm.rank() == 0 && verbose == true){
        mu.setIcntl(1, 6);
        mu.setIcntl(2, 6);
        mu.setIcntl(3, 6);
    }
    dmumps_c(&mu);
    t = MPI_Wtime() - t;

    if(inter_comm.rank() == 0){
        mu.setIcntl(1, -1);
        mu.setIcntl(2, -1);
        mu.setIcntl(3, -1);
    }

    if(mu.getInfo(1) != 0){
        cout << string(32, '-') << endl
             << "| MUMPS Factoriz FAILED on MA " << setw(7) << inter_comm.rank() << " |" << endl
             << string(32, '-') << endl
             << "| info(1)       : " << setw(6) << mu.getInfo(1) << string(4, ' ') << " |" << endl
             << "| info(2)       : " << setw(6) << mu.getInfo(2) << string(4, ' ') << " |" << endl
             << string(32, '-') << endl;
        throw 100 * mu.getInfo(1) - mu.getInfo(2);
    }

    double smem;
    if(instance_type == 0) {
        double mem = mu.getInfoG(21);

        if(IRANK == 0) mpi::reduce(inter_comm, mem, smem, std::plus<double>(), 0);
        else mpi::reduce(inter_comm, mem, std::plus<double>(), 0);

        if(IRANK == 0) smem = smem/world.size();
    }

    if(instance_type == 0 && verbose == true) {
        double flop = mu.getRinfoG(3);
        int prec = cout.precision();
        cout.precision(2);
        cout << string(32, '-') << endl
             << "| MUMPS FACTORIZ on MA " << setw(7) << inter_comm.rank() << " |" << endl
             << string(32, '-') << endl
             << "| N             : " << setw(12) << mumps.n << " |" << endl
             << "| NZ            : " << setw(12) << mumps.nz << " |" << endl
             << "| Flops         : " << setw(6) << scientific << flop << string(4, ' ') << " |" << endl
             << "| Time          : " << setw(6) << t << " sec |" << endl
             << "| avg memory    : " << setw(6) << smem << " M| "<< endl
             << string(32, '-') << endl;
        cout.precision(prec);
    }
}
