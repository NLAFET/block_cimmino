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
#include "paircomparison.h"

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

        // choice 0: slave assigned as encountered
        // choice 1:
        //    a) for each master in descending order of flops: put slaves in same node
        //    b) for each master in descending order of flops: add remaining slaves GROUPED in remaining slots
        // choice 2: idem with merged loops
        //    a) for each master in descending order of flops: fill node then add remaining
        // choice -1: previous implementation (masters first ranks: no longer working)
        int choice=icntl[Controls::slave_def];
        //get correspondance between comm and inter_comm ranks for masters
        std::vector<int> masters_comm_rank;
        if (choice >= 0) {
            mpi::gather(inter_comm, comm.rank(), masters_comm_rank, 0);
        }
        // Sort masters per #slaves via tmp_ord to ord_nodes
        std::vector<int> ord_masters; // masters ordered index vs. original index
        if (choice > 0) {
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
        }
        // DEFINE SLAVES
        if (choice==0) {
            // distribute slaves to masters according to previous count
            int current_master=0;
            for(int proc=0; proc<instance_type_vect.size(); ++proc) {
                if (instance_type_vect[proc]) {
                    // let root MPI inform slaves of their master
                    if(comm.rank() == 0) {
                        // current_master is numbered as in inter_comm but should be sent as in comm
                        // with the previous implementation there was no difference
                        comm.send(proc, 11, masters_comm_rank[current_master]);
                    }
                    // masters recognize their slaves
                    if(inter_comm.rank() == current_master) {
                        my_slaves.push_back(proc);
                    }
                    // when a master has all its slaves, go to next
                    --slaves_for_me[current_master];
                    if (!slaves_for_me[current_master])
                        ++current_master;
                }
            }
        } else if (choice==1) { // first choice
            if (comm.rank() == 0) {
                for(int idx = 0 ; idx < ord_masters.size(); idx++) {
                    std::vector<int> tmp_my_slaves;
                    int master=ord_masters[idx];
//                    LINFO << "MASTER: " << master;
                    int node_master=masters_node[master];
                    int node=node_master;
                    while (slaves_for_me[master] > 0 && node_map_slaves[node_master].size() > 0) {
                        int slave=node_map_slaves[node][0];
                        comm.send(slave, 11, masters_comm_rank[master]);
                        tmp_my_slaves.push_back(slave);
                        node_map_slaves[node].erase(node_map_slaves[node].begin());
                        --slaves_for_me[master];
                    }
                    if (masters_comm_rank[master] == 0) {
                        my_slaves=tmp_my_slaves;
                    } else {
                        inter_comm.send(master, 63, tmp_my_slaves);
                    }
                }
//                LINFO << "SECOND STEP";
                for(int idx = 0 ; idx < ord_masters.size(); idx++) {
                    std::vector<int> tmp_my_slaves;
                    int master=ord_masters[idx];
//                    LINFO << "MASTER: " << master;
                    int node_master=masters_node[master];
                    int iii=0;
                    int jjj=0;
                    while (slaves_for_me[master] > 0) {
                        int slave=node_map_slaves[iii][0];
                        comm.send(slave, 11, masters_comm_rank[master]);
                        tmp_my_slaves.push_back(slave);
                        node_map_slaves[iii].erase(node_map_slaves[iii].begin());
                        --slaves_for_me[master];

                        ++jjj;
                        if (jjj==node_map_slaves[iii].size()) {
                            jjj=0;
                            ++iii;
                        }
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
        } else if (choice==2) { // second choice
            if (comm.rank() == 0) {
                for(int idx = 0 ; idx < ord_masters.size(); idx++) {
                    std::vector<int> tmp_my_slaves;
                    int master=ord_masters[idx];
//                    LINFO << "MASTER: " << master;
                    int node_master=masters_node[master];
                    int node=node_master;
                    int iii=0;
                    int jjj=0;
                    while (slaves_for_me[master] > 0) {
                        if (node_map_slaves[node_master].size() == 0) node=iii;
                        int slave=node_map_slaves[node][0];
                        comm.send(slave, 11, masters_comm_rank[master]);
                        tmp_my_slaves.push_back(slave);
                        node_map_slaves[node].erase(node_map_slaves[node].begin());
                        --slaves_for_me[master];
                        if (node_map_slaves[node].size() == 0) {
                            ++jjj;
                            if (jjj==node_map_slaves[node].size()) {
                                jjj=0;
                                ++iii;
                            }
                        }
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
        } else {
            int current_slave = inter_comm.size();
            for(int your_master = 0 ; your_master < inter_comm.size(); your_master++) {
                for(int i = 0; i < slaves_for_me[your_master]; i++) {
                    // let root MPI inform slaves of their master
                    if(inter_comm.rank() == 0)
                        comm.send(current_slave, 11,  your_master);

                    // masters recognize their slaves
                    if(inter_comm.rank() == your_master) {
                        my_slaves.push_back(current_slave);
                    }
                    current_slave++;
                }
            }
        }

/*        std::string res = "";
        res.append("MASTER ");
        res.append(std::to_string(comm.rank()));
        res.append(": SLAVES ");
        for(int jjj=0; jjj<my_slaves.size(); ++jjj) {
            res.append(std::to_string(my_slaves[jjj]));
            res.append("; ");
        }
        std::cout << res << "\n";*/
        // Now that the slaves know who's their daddy, daddy tells who are their brothers
        for(std::vector<int>::iterator slave = my_slaves.begin(); slave != my_slaves.end(); ++slave) {
            comm.send(*slave, 12, my_slaves);
        }
    } else {
        comm.recv(0, 11, my_master);
        //Store in my_slaves my brothers in slavery
        comm.recv(my_master, 12, my_slaves);
    }

    LINFO2 << "Process " << inter_comm.rank()
           << " has " << my_slaves.size()
           << " workers";
}
