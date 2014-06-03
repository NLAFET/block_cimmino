#include <abcd.h>
#include <mumps.h>

void abcd::analyseAugmentedSystems(MUMPS &mu)
{

  double t = MPI_Wtime();
  mu(1);

  t = MPI_Wtime() - t;

  if(mu.getInfo(1) != 0){
    LERROR << string(32, '-') ;
    LERROR << "| MUMPS ANALYSIS FAILED on MA " << setw(7) << inter_comm.rank() << " |" ;
    LERROR << string(32, '-') ;
    LERROR << "| info(1)       : " << setw(6) << mu.getInfo(1) << string(4, ' ') << " |" ;
    LERROR << "| info(2)       : " << setw(6) << mu.getInfo(2) << string(4, ' ') << " |" ;
    LERROR << string(32, '-') ;;

    LERROR << "MUMPS exited with " << mumps_S.info[0];
    int job = -70 + mu.info[0];
    mpi::broadcast(intra_comm, job, 0);
    throw std::runtime_error("MUMPS exited with an error");
  }

  if(instance_type == 0) {
    double flop = mu.getRinfo(1);
    int prec = cout.precision();
    cout.precision(2);
    LINFO << string(32, '-') ;
    LINFO << "| MUMPS ANALYSIS on MA " << setw(7) << inter_comm.rank() << " |" ;
    LINFO << string(32, '-') ;
    LINFO << "| Flops estimate: " << setw(6) << scientific << flop << string(4, ' ') << " |" ;
    LINFO << "| Time:           " << setw(6) << t << " sec |" ;
    LINFO << string(32, '-') ;;
    cout.precision(prec);
  }
}
