#include<abcd.h>
#include<mumps.h>

void abcd::factorizeAugmentedSystems(MUMPS &mu)
{
    mpi::communicator comm;
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
        LERROR << string(32, '-') ;
        LERROR << "> MUMPS Factoriz FAILED on master " << inter_comm.rank();
        LERROR << "MUMPS exited with " << mu.getInfo(1);
        
        int job = -700 + mu.getInfo(1);
        mpi::broadcast(intra_comm, job, 0);

        // oh the infamous -9!
        if (mu.getInfo(1) == -9) {
            LERROR << "MUMPS's internet real workarray is too small.";
            LERROR << "MUMPS is missing "
                   << setprecision(2)
                   << (mu.getInfo(2) > 0 ? mu.getInfo(2) * sizeof(double) : mu.getInfo(2) * pow(10, 6) * sizeof(double))
                   << " Bytes";
            LERROR << "MUMPS ICNTL[14] = " << mu.getIcntl(14);

        } else if (mu.getInfo(1) == -10) {
            LERROR << "MUMPS says that the augmented matrix is singular.";
            LERROR << "That should not happen if your system is full column rank";
        }
        
        throw std::runtime_error("MUMPS exited with an error!");
    }

    double smem;
    if(instance_type == 0) {
        double mem = mu.getInfoG(21);

        if(IRANK == 0) mpi::reduce(inter_comm, mem, smem, std::plus<double>(), 0);
        else mpi::reduce(inter_comm, mem, std::plus<double>(), 0);

        if(IRANK == 0) smem = smem/inter_comm.size();
    }

    if(instance_type == 0) {
        double flop = mu.getRinfoG(3);
        int prec = cout.precision();
        cout.precision(2);
        LINFO << string(32, '-') ;
        LINFO << "| MUMPS FACTORIZ on MA " << setw(7) << inter_comm.rank() << " |" ;
        LINFO << string(32, '-') ;
        LINFO << "| N             : " << setw(12) << mumps.n << " |" ;
        LINFO << "| NZ            : " << setw(12) << mumps.nz << " |" ;
        LINFO << "| Flops         : " << setw(6) << scientific << flop << string(4, ' ') << " |" ;
        LINFO << "| Time          : " << setw(6) << t << " sec |" ;
        LINFO << "| avg memory    : " << setw(6) << smem << " M| ";
        LINFO << string(32, '-') ;;
        cout.precision(prec);
    }
}
