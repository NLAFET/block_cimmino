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
        LERROR << "MUMPS exited with TWO " << mu.getInfo(2);
        
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
    int max_mem;
    int sflop;
    if(instance_type == 0) {
        double mem = mu.getInfoG(22)/intra_comm.size();
        mpi::reduce(inter_comm, mem, smem, std::plus<double>(), 0);
        if(IRANK == 0) smem = smem/inter_comm.size();

        int mmem = mu.getInfoG(21);
	mpi::reduce(inter_comm, mmem, max_mem, mpi::maximum<int>(),0);

        int flop = mu.getRinfoG(3);
	mpi::reduce(inter_comm, flop, sflop, std::plus<int>(),0);

        int prec = cout.precision();
        cout.precision(2);
        LINFO << string(32, '-') ;
        LINFO << "| MUMPS FACTORIZ on MA " << setw(7) << inter_comm.rank() << " |" ;
        LINFO << string(32, '-') ;
        LINFO << "| N             : " << setw(12) << mumps.n << " |" ;
        LINFO << "| NZ            : " << setw(12) << mumps.nz << " |" ;
        LINFO << "| Flops         : " << setw(6) << scientific << flop << string(4, ' ') << " |" ;
        LINFO << "| Time          : " << setw(6) << t << " sec |" ;
        LINFO << "| Average memory    : " << setw(6) << mem << " M| ";
        LINFO << string(32, '-') ;;
        cout.precision(prec);
    }

    if(IRANK == 0) LINFO << "Factorization average memory : " << setw(6) << smem << " M";
    if(IRANK == 0) LINFO << "Factorization maximum memory : " << setw(6) << max_mem << " M";
    if(IRANK == 0) LINFO << "Factorization total flops : " << setw(6) << sflop << " flops";
}
