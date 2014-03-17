#include <abcd.h>
#include <fstream>

void abcd::initializeCimmino()
{
    mpi::communicator world;
    int *sym_perm;

    mpi::broadcast(world, nbparts, 0);

    
    if(world.size() > parallel_cg) {
        if(instance_type == 0) {
            if(inter_comm.rank() == 0 && instance_type == 0)
                cout << "[+] Initializing MUMPS" << endl;
            initializeMumps(mumps, true);
            createAugmentedSystems(n_aug, nz_aug, irn_aug, jcn_aug, val_aug);

            mumps.n = n_aug;
            mumps.nz = nz_aug;
            mumps.irn = &irn_aug[0];
            mumps.jcn = &jcn_aug[0];
            mumps.a = &val_aug[0];

            if(inter_comm.rank() == 0 && instance_type == 0)
                cout << "[+] Launching Initial MUMPS analysis" << endl;
            analyseAugmentedSystems(mumps);
            
            sym_perm = new int[mumps.n];
            std::copy(mumps.sym_perm, mumps.sym_perm + mumps.n, sym_perm);

        }


        if(instance_type == 0) {
            mumps.job = -2;
            dmumps_c(&mumps);
        }

        allocateMumpsSlaves(mumps);
        initializeMumps(mumps);

    } else {
        allocateMumpsSlaves(mumps);
        initializeMumps(mumps);

        createAugmentedSystems(n_aug, nz_aug, irn_aug, jcn_aug, val_aug);

    }

    if(instance_type == 0) {
        mumps.n = n_aug;
        mumps.nz = nz_aug;
        mumps.irn = &irn_aug[0];
        mumps.jcn = &jcn_aug[0];
        mumps.a = &val_aug[0];
    }
    
    if(inter_comm.rank() == 0 && instance_type == 0)
        cout << "[+] Launching MUMPS analysis" << endl;
    abcd::analyseAugmentedSystems(mumps);

    //delete[] sym_perm;

}
