#include <abcd.h>
#include <mumps.h>

void abcd::analyseAugmentedSystems(MUMPS &mu)
{
    mu.job = 1;

    double t = MPI_Wtime()

    dmumps_c(&mu);

    t = MPI_Wtime() - t;

    if(mu.getInfo(1) != 0){
        cout << string(32, '-') << endl
             << "| MUMPS ANALYSIS FAILED on MA " << setw(7) << inter_comm.rank() << " |" << endl
             << string(32, '-') << endl
             << "| info(1)       : " << setw(6) << mu.getInfo(1) << string(4, ' ') << " |" << endl
             << "| info(2)       : " << setw(6) << mu.getInfo(2) << string(4, ' ') << " |" << endl
             << string(32, '-') << endl;
        throw 100 * mu.getInfo(1) - mu.getInfo(2);
    }

    if(instance_type == 0 && verbose == true) {
        double flop = mu.getRinfo(1);
        int prec = cout.precision();
        cout.precision(2);
        cout << string(32, '-') << endl
             << "| MUMPS ANALYSIS on MA " << setw(7) << inter_comm.rank() << " |" << endl
             << string(32, '-') << endl
             << "| Flops estimate: " << setw(6) << scientific << flop << string(4, ' ') << " |" << endl
             << "| Time:           " << setw(6) << t << " sec |" << endl
             << string(32, '-') << endl;
        cout.precision(prec);
    }
}
