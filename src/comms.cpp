#include "abcd.h"

/// Assignes each mpi-process to its category : CG-master or MUMPS-Slave
void abcd::inter_group_mapping()
{
    mpi::communicator world;
    mpi::group grp = world.group();

    mpi::communicator cm;
    mpi::group gp;

    mpi::broadcast(world, parallel_cg, 0);

    if(parallel_cg > world.size()) throw -14;

    inter_comm = world.split(world.rank() < parallel_cg);
    
    if (world.rank() < parallel_cg )
        instance_type = 0;
    else
        instance_type = 1;
}

/// Distribute the partitions over CG processes
void abcd::distribute_parts()
{
    
}