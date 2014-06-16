#include <abcd.h>
#include <fstream>

void abcd::initializeCimmino()
{
    int *sym_perm;

    mpi::broadcast(comm, icntl[Controls::nbparts], 0);
    
    if(comm.size() > parallel_cg) {
        if(instance_type == 0) {
            if(inter_comm.rank() == 0 && instance_type == 0)
                LINFO << "Initializing MUMPS";
            initializeMumps(mumps, true);
            createAugmentedSystems(n_aug, nz_aug, irn_aug, jcn_aug, val_aug);

            mumps.n = n_aug;
            mumps.nz = nz_aug;
            mumps.irn = &irn_aug[0];
            mumps.jcn = &jcn_aug[0];
            mumps.a = &val_aug[0];

            if(inter_comm.rank() == 0 && instance_type == 0)
                LINFO << "Launching Initial MUMPS analysis";
            analyseAugmentedSystems(mumps);
            
            sym_perm = new int[mumps.n];
            std::copy(mumps.sym_perm, mumps.sym_perm + mumps.n, sym_perm);

        }


        if(instance_type == 0) {
            mumps.job = -2;
            dmumps_c(&mumps);

            mumps.initialized = false;
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
}
