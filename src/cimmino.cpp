#include "abcd.h"
#include <fstream>

void abcd::initializeCimmino()
{
    mpi::communicator world;
    if(instance_type == 0) {
        initializeMumps(true);
        createAugmentedSystems();
        if(world.rank() == 0) {
            ofstream myfile;
            myfile.open("/tmp/m");
            myfile << "%%MatrixMarket matrix coordinate real general" << endl;
            myfile << mumps.n << " " << mumps.n << " " << mumps.nz << endl;

            for(int i = 0; i < mumps.nz; i++) {
                myfile << mumps.irn[i] << " " << mumps.jcn[i] << " " << scientific << mumps.a[i] << endl;
            }
            myfile.close();
        }
        cout << "[+] Launching MUMPS analysis" << endl;
        analyseAugmentedSystems();
    }

    allocateMumpsSlaves();

    if(instance_type == 0) {
        mumps.job = -2;
        dmumps_c(&mumps);
    }

    world.barrier();
    initializeMumps();

    if(instance_type == 0) {
        createAugmentedSystems();
    }

}
