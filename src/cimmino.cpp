#include <abcd.h>
#include <fstream>

void abcd::initializeCimmino()
{
    mpi::communicator world;
    int *sym_perm;
    
    if(instance_type == 0) {
        if(inter_comm.rank() == 0 && instance_type == 0)
            cout << "[+] Initializing MUMPS" << endl;
        initializeMumps(true);
        createAugmentedSystems();
/*         if(world.rank() == 0) {
 *             ofstream myfile;
 *             myfile.open("/tmp/m");
 *             myfile << "%%MatrixMarket matrix coordinate real general" << endl;
 *             myfile << mumps.n << " " << mumps.n << " " << mumps.nz << endl;
 * 
 *             for(int i = 0; i < mumps.nz; i++) {
 *                 myfile << mumps.irn[i] << " " << mumps.jcn[i] << " " << scientific << mumps.a[i] << endl;
 *             }
 *             myfile.close();
 *         }
 */
        if(inter_comm.rank() == 0 && instance_type == 0)
            cout << "[+] Launching Initial MUMPS analysis" << endl;
        analyseAugmentedSystems();
        
        sym_perm = new int[mumps.n];
        std::copy(mumps.sym_perm, mumps.sym_perm + mumps.n, sym_perm);

    }

    allocateMumpsSlaves();

    if(instance_type == 0) {
        mumps.job = -2;
        dmumps_c(&mumps);
    }

    world.barrier();
    initializeMumps();
    //ordering given
    //setMumpsIcntl(7,1);
    //setMumpsIcntl(6,5);
    //setMumpsIcntl(8,-2);
    //setMumpsIcntl(28,2);

    if(instance_type == 0) {
        createAugmentedSystems();
        //mumps.perm_in = sym_perm;
    }
    
    if(inter_comm.rank() == 0 && instance_type == 0)
        cout << "[+] Launching MUMPS analysis" << endl;
    abcd::analyseAugmentedSystems();

}
