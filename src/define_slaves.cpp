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

#include <abcd.h>

void abcd::allocateMumpsSlaves(MUMPS &mu)
{
    if(instance_type == 0) {
        std::vector<long> flops(inter_comm.size()); // MUMPS estimated flops for each MPI
        mpi::all_gather(inter_comm, (long) mu.getRinfo(1), flops);
        std::vector<int> slaves_for_me(inter_comm.size()); // number of slaves per master

        // each slave is assigned to the master with current max flops
        // each time a slave is assigned, its master's flops becomes the original flops divided by #slaves+1
	int nb_slaves = comm.size() - inter_comm.size();
	std::vector<long> flops_orig(flops); // original vector of flops
	for(int i = 0; i < nb_slaves; i++){
	    int max_index = max_element_index(flops.begin(), flops.end());
	    slaves_for_me[max_index]++;
            flops[max_index] = flops_orig[max_index] / (slaves_for_me[max_index] + 1);
	}

        // Sort masters per #slaves via tmp_ord to ord_nodes
        std::vector<int> ord_masters; // masters ordered index vs. original index
        std::vector<std::pair<int, int>> tmp_ord; // tmp vector to sort masters by #slaves
        for (int iii=0; iii<slaves_for_me.size(); ++iii) {
            std::pair<int, int> p (iii, slaves_for_me[iii]);
            tmp_ord.push_back(p);
        }
        std::sort(tmp_ord.begin(), tmp_ord.end(),
            pair_comparison<int, int,
                position_in_pair::second,
                comparison_direction::descending>);
        for (int iii=0; iii<tmp_ord.size(); ++iii) {
            ord_masters.push_back(tmp_ord[iii].first);
        }

        //get correspondance between comm and inter_comm ranks for masters
        std::vector<int> masters_comm_rank;
        std::vector<std::vector<int>> disp_node_map_slaves; // to display info about slaves
        std::vector<std::vector<int>> disp_node_map_masters; // to display info about masters
        mpi::gather(inter_comm, comm.rank(), masters_comm_rank, 0);
        if (!comm.rank()) {
            for (int i=0; i<node_map_slaves.size(); ++i) {
                std::vector<int> v;
                disp_node_map_slaves.push_back(v);
                disp_node_map_masters.push_back(v);
            }
            for (int i=0; i<ord_masters.size(); ++i) {
                int master=ord_masters[i];
                disp_node_map_masters[masters_node[master]].push_back(masters_comm_rank[master]);
            }
        }

        int choice=icntl[Controls::slave_def];
        // DEFINE SLAVES
        // choice 0: slave assigned as encountered
        if (choice==0) {
            // distribute slaves to masters according to previous count
            int current_master=0;
            for(int proc=0; proc<instance_type_vect.size(); ++proc) {
                if (instance_type_vect[proc]) {
                    // let root MPI inform slaves of their master
                    if(!comm.rank()) {
                        // current_master is numbered as in inter_comm but should be sent as in comm
                        // with the previous implementation there was no difference
                        comm.send(proc, 11, masters_comm_rank[current_master]);
                    }
                    // masters recognize their slaves
                    if(inter_comm.rank() == current_master)
                        my_slaves.push_back(proc);
                    // when a master has all its slaves, go to next
                    --slaves_for_me[current_master];
                    if(!comm.rank())
                        disp_node_map_slaves[mpi_map[proc]].push_back(masters_comm_rank[current_master]);
                    if (!slaves_for_me[current_master])
                        ++current_master;
                }
            }
        // choice 1:
        //    a) for each master in descending order of flops: put slaves in same node
        //    b) for each master in descending order of flops: add remaining slaves GROUPED in remaining slots
        } else if (choice==1) { // first choice
            if (comm.rank() == 0) {
                for(int idx = 0 ; idx < ord_masters.size(); idx++) {
                    std::vector<int> tmp_my_slaves;
                    int master=ord_masters[idx];
                    int node_master=masters_node[master];
                    int node=node_master;
                    while (slaves_for_me[master] > 0 && node_map_slaves[node_master].size() > 0) {
                        int slave=node_map_slaves[node][0];
                        comm.send(slave, 11, masters_comm_rank[master]);
                        tmp_my_slaves.push_back(slave);
                        node_map_slaves[node].erase(node_map_slaves[node].begin());
                        --slaves_for_me[master];
                        disp_node_map_slaves[node].push_back(masters_comm_rank[master]);
                    }
                    if (masters_comm_rank[master] == 0) {
                        my_slaves=tmp_my_slaves;
                    } else {
                        inter_comm.send(master, 63, tmp_my_slaves);
                    }
                }
                for(int idx = 0 ; idx < ord_masters.size(); idx++) {
                    std::vector<int> tmp_my_slaves;
                    int master=ord_masters[idx];
                    int node_master=masters_node[master];
                    int iii=0;
                    while (slaves_for_me[master] > 0) {
                        // Go to the next non-empty node if the current one is empty
                        while(node_map_slaves[iii].size() == 0)
                           ++iii;
                        int slave=node_map_slaves[iii][0];
                        comm.send(slave, 11, masters_comm_rank[master]);
                        tmp_my_slaves.push_back(slave);
                        node_map_slaves[iii].erase(node_map_slaves[iii].begin());
                        --slaves_for_me[master];
                        disp_node_map_slaves[iii].push_back(masters_comm_rank[master]);
                    }
                    if (masters_comm_rank[master] == 0) {
                        my_slaves.insert(my_slaves.end(), tmp_my_slaves.begin(), tmp_my_slaves.end());
                    } else {
                        inter_comm.send(master, 64, tmp_my_slaves);
                    }
                }
            } else {
                inter_comm.recv(0, 63, my_slaves);
                std::vector<int> tmp_my_slaves;
                inter_comm.recv(0, 64, tmp_my_slaves);
                my_slaves.insert(my_slaves.end(), tmp_my_slaves.begin(), tmp_my_slaves.end());
            }
        // choice 2: idem with merged loops
        //    a) for each master in descending order of flops: fill node then add remaining
        } else if (choice==2) { // second choice
            if (comm.rank() == 0) {
                for(int idx = 0 ; idx < ord_masters.size(); idx++) {
                    std::vector<int> tmp_my_slaves;
                    int master=ord_masters[idx];
                    int node_master=masters_node[master];
                    int node=node_master;
                    int iii=0;
                    while (slaves_for_me[master] > 0) {
                        // If the node of the master and the current node are empty,
                        // go to the next non-empty one
                        if (node_map_slaves[node_master].size() == 0) {
                            while(node_map_slaves[iii].size() == 0)
                                ++iii;
                            node=iii;
                        }
                        int slave=node_map_slaves[node][0];
                        comm.send(slave, 11, masters_comm_rank[master]);
                        tmp_my_slaves.push_back(slave);
                        node_map_slaves[node].erase(node_map_slaves[node].begin());
                        --slaves_for_me[master];
                        disp_node_map_slaves[node].push_back(masters_comm_rank[master]);
                    }
                    if (masters_comm_rank[master] == 0) {
                        my_slaves=tmp_my_slaves;
                    } else {
                        inter_comm.send(master, 63, tmp_my_slaves);
                    }
                }
            } else {
                inter_comm.recv(0, 63, my_slaves);
            }
        // choice -1: previous implementation (masters first ranks: no longer working)
        } else {
            int current_slave = inter_comm.size();
            for(int your_master = 0 ; your_master < inter_comm.size(); your_master++) {
                for(int i = 0; i < slaves_for_me[your_master]; i++) {
                    // let root MPI inform slaves of their master
                    if(inter_comm.rank() == 0)
                        comm.send(current_slave, 11,  your_master);
                    if(!comm.rank())
                        disp_node_map_slaves[mpi_map[current_slave]].push_back(masters_comm_rank[your_master]);

                    // masters recognize their slaves
                    if(inter_comm.rank() == your_master)
                        my_slaves.push_back(current_slave);
                    current_slave++;
                }
            }
        }

        // Now that the slaves know who's their daddy, daddy tells who are their brothers
        for(std::vector<int>::iterator slave = my_slaves.begin(); slave != my_slaves.end(); ++slave) {
            comm.send(*slave, 12, my_slaves);
        }

        // Display Nodes/Masters/Slaves information
        if (!comm.rank()) {
            stringstream res;
            for(int node=0; node<disp_node_map_masters.size(); ++node) {
                res << "Node " << node << ": Masters ";
                for(int master=0; master<disp_node_map_masters[node].size(); ++master) {
                    res << disp_node_map_masters[node][master];
                    if(master!=disp_node_map_masters[node].size()-1)
                        res << ",";
                }
                res << " with slaves of masters: ";
                for(int slave=0; slave<disp_node_map_slaves[node].size(); ++slave) {
                    res << disp_node_map_slaves[node][slave];
                    if(slave!=disp_node_map_slaves[node].size()-1)
                        res << ",";
                }
                LINFO2 << res.str();
                res.str("");
            }
        }
    } else {
        comm.recv(0, 11, my_master);
        //Store in my_slaves my brothers in slavery
        comm.recv(my_master, 12, my_slaves);
    }
//    LINFO2 << "Process " << inter_comm.rank()
    LINFO2 << "Process " << comm.rank()
        << " has " << my_slaves.size()
        << " workers";
}
