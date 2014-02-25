#include "abcd.h"
#include "mumps.h"

#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/exceptions.hpp>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>


using namespace std;
using namespace boost::property_tree;

int main(int argc, char* argv[]) 
{
    mpi::environment env(argc, argv);
    mpi::communicator world;

    abcd obj;

    // This should be done only by the master
    if(world.rank() == 0) {
        string config_file;
        if(argc < 2) {
            config_file = "config_file.info";
            //clog << "Usage " << argv[0] << " config_file.json" << endl;
            //exit(-1);
        } else if(argc >= 2){
            config_file = argv[1];
        }

        /* open the config file */
        ifstream conf_file(config_file.c_str());
        boost::property_tree::ptree pt;

        if(conf_file){
            try {
                read_info(conf_file, pt);
            } catch (info_parser::info_parser_error e){
                clog << " *** Problem while parsing the info file" << endl
                     << " *** The what() : " << e.what() << endl << endl
                     << " *   Check that there is no aditional comma at the end" << endl
                     << " *   or an unamed node" << endl;
                exit(-1);
            }
            conf_file.close();
        } else {
            clog << "Error while opening config file" << endl;
            exit(-2);
        }

        /* READING THE MATRIX AND THE RHS */
        string matrix_file;
        try {
            if(argc <= 2) matrix_file = pt.get<string>("system.matrix_file");
            else matrix_file = argv[2];
        } catch (ptree_bad_path e) {
            clog << "Error parsing the file, you have to give the matrix name as in the example file" << endl;
            clog << "The what() : " << e.what() << endl << endl;
            exit(-1);
        }

        FILE *f = fopen(matrix_file.c_str(), "r");

        if(f == NULL){
            cerr << "Error opening the file '"<< matrix_file << "'" << endl;
            exit(-1);
        }
        MM_typecode mat_code;

        mm_read_banner(f, &mat_code);
        mm_read_mtx_crd_size(f, (int *)&obj.m, (int *)&obj.n, (int *)&obj.nz);

        if(mm_is_symmetric(mat_code))
            obj.sym = true;
        else
            obj.sym = false;

        // allocate the arrays
        obj.irn = new int[obj.nz];
        obj.jcn = new int[obj.nz];
        obj.val = new double[obj.nz];

        mm_read_mtx_crd_data(f, obj.m, obj.n, obj.nz, obj.irn, obj.jcn, obj.val, mat_code);

        cout << "Matrix information : ";
        cout << "m = " << obj.m << "; n = " << obj.n << "; nz = " << obj.nz << endl;

        fclose(f);

        //read the rhs here!
        boost::optional<string> rhs_file = pt.get_optional<string>("system.rhs_file");
        
        if(rhs_file){
            FILE *rhs_f = fopen(rhs_file->c_str(), "r");

            if(rhs_f == NULL){
                cerr << "Error opening the file '"<< *rhs_file << "'" << endl;
                exit(-1);
            }

            MM_typecode rhs_code;
            int nb_v, m_v;

            mm_read_banner(rhs_f, &rhs_code);
            mm_read_mtx_array_size(rhs_f, &m_v, &nb_v);

            obj.rhs = new double[nb_v * m_v];
            cout << "Reading "<< m_v << " values for each of the " << nb_v << " rhs" << endl;

            double cv;
            for(int i = 0; i < m_v*nb_v; i++){
                fscanf(rhs_f, "%lf", &cv);
                obj.rhs[i] = cv;
            }
            //
            fclose(rhs_f);
        }

        bool testMumps = pt.get<bool>("testMumps", false);
        MUMPS mu;
        if(testMumps){

            if(obj.sym) {
                mu.sym = 2;
            } else {
                mu.sym = 0;
            }
            mu.par = 1;
            mu.job = -1;

            if (obj.rhs != NULL){
                int todo = 6;
                mpi::broadcast(world, todo, 0);
            } else {
                int todo = 4;
                mpi::broadcast(world, todo, 0);
            }
            mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) world);

            dmumps_c(&mu);

            mu.setIcntl(1, -1);
            mu.setIcntl(2, -1);
            mu.setIcntl(3, -1);

            mu.setIcntl(6, 5);
            mu.setIcntl(7, 5);
            mu.setIcntl(8, -2);
            mu.setIcntl(12, 2);
            mu.setIcntl(14, 90);

            mu.n = obj.n;
            mu.nz = obj.nz;
            mu.irn = obj.irn;
            mu.jcn = obj.jcn;
            mu.a = obj.val;

            if (obj.rhs != NULL){
                mu.job = 6;
                mu.nrhs = 1;
                mu.rhs = new double[mu.n];
                for(int i = 0; i < mu.n; i++)
                    mu.rhs[i] = obj.rhs[i];
            } else {
                mu.job = 4;
            }

            cout << "=============================" << endl;
            cout << "LAUNCHING MUMPS ON THE MATRIX" << endl;
            cout << "-----------------------------" << endl;
            double t = MPI_Wtime();
            dmumps_c(&mu);
            clog << "MUMPS RUN FINISHED in " << MPI_Wtime() - t << endl;
            clog << "=============================" << endl;
        } else {
            testMumps = 0;
            mpi::broadcast(world, testMumps, 0);
        }

        //
        /*  DONE READING THE MATRIX AND THE RHS */

        try{
            ptree::key_type partitioning = pt.get<ptree::key_type>("partitioning");
        } catch (ptree_bad_path e){
            clog << "Error parsing the file, you have to give the partitioning as in the example file" << endl;
            clog << "The what() : " << e.what() << endl << endl;
            exit(-1);
        }

        if(argc <= 3) obj.nbparts = pt.get<int>("partitioning.nbparts");
        else obj.nbparts = atoi(argv[3]);

        if(obj.nbparts < 0 && obj.nbparts <= obj.m) {
            clog << "Error parsing the file, the number of partitions has to be positive and smaller than the number of rows" << endl;
            exit(-1);
        }

        obj.partitioning_type = pt.get<int>("partitioning.type", 2);
        obj.guessPartitionsNumber = pt.get<int>("partitioning.guess", 0);
        obj.dcntl[8] = pt.get<double>("partitioning.imba", 0.5);

        if(obj.partitioning_type == 1){
            string parts = pt.get<string>("partitioning.partsfile", "");
            std::vector<int> nrows;

            if(parts.length() == 0){
                try{
                    ptree::key_type partitioning = pt.get<ptree::key_type>("partitioning.nbrows");
                } catch (ptree_bad_path e){
                    clog << "Error parsing the file, you have to give the number of rows per partition as in the example file" << endl;
                    clog << "The what() : " << e.what() << endl << endl;
                    exit(-1);
                }

                BOOST_FOREACH( ptree::value_type v, pt.get_child("partitioning.nbrows") )
                {
                    nrows.push_back(atoi(v.first.data()));
                }

            } else {
                ifstream f;
                f.open(parts.c_str());

                for(unsigned k = 0; k < (unsigned int)obj.nbparts; k++) {
                    string l;
                    getline(f, l);
                    nrows.push_back(atoi(l.c_str()));
                }

                f.close();
            }

            if(nrows.size() != (size_t)obj.nbparts){
                clog << "Error parsing the file, nbparts is different from the partitioning description" << endl;
                exit(-1);
            }

            obj.nbrows = VECTOR_int(&nrows[0], obj.nbparts);
        }

        obj.write_problem   = pt.get<string>("write_problem", "");

        obj.icntl[8]    = pt.get<int>("esparse", 0);
        obj.icntl[9]    = pt.get<int>("scaling", 2);

        boost::optional<ptree::key_type> augmentation = pt.get_optional<ptree::key_type>("augmentation");

        if(augmentation){
            obj.icntl[10]   = pt.get<int>("augmentation.type", 2);
            obj.dcntl[10]   = pt.get<double>("augmentation.filtering", 0.0);
            obj.icntl[11]   = pt.get<int>("augmentation.analysis", 0);
            obj.icntl[12]   = pt.get<int>("augmentation.project_only", 0);
            obj.icntl[13]   = pt.get<int>("augmentation.denserhs", 0);
            obj.icntl[14]   = pt.get<int>("augmentation.multirhs", 256);
            obj.icntl[15]   = pt.get<int>("augmentation.iterative", 0);
            obj.dcntl[15]   = pt.get<double>("augmentation.precond", 0.0);
            obj.write_s     = pt.get<string>("augmentation.write_s", "");
        }

        obj.parallel_cg = pt.get<int>("dist_scheme", obj.nbparts < world.size() ? obj.nbparts : world.size());

        bool error = false;
        try {
            

            obj.verbose =  pt.get<int>("all_verbose", 0);

            obj.bc(-1);

            double t = MPI_Wtime();
            obj.bc(1);
            obj.bc(2);

            obj.nrhs = 1;
            obj.use_xf = false;

            if(argc <= 4) obj.block_size = pt.get<int>("system.block_size", 1);
            else obj.block_size = atoi(argv[4]);

            obj.itmax = pt.get<int>("system.itmax", 2000);
            obj.threshold = pt.get<double>("system.threshold", 1e-12);

            obj.verbose =  pt.get<int>("solve_verbose", 0);
            obj.bc(3);

            clog << "Total time: " << MPI_Wtime() - t << endl;
        } catch(int e) {
            cout << world.rank() << " Error code : " << e << endl;
            mpi::broadcast(world, e, 0);
            error = true;
        }

        if(testMumps && !error) {
            double infTop = 0, infBot = 0;
            if(mu.job == 6){
                ofstream f; 
                f.open("/tmp/out_comp");
                for(int i = 0; i < mu.n; i++) {
                    f << mu.rhs[i] << "\t" << obj.sol(i,0) << "\n";
                    if (abs(mu.rhs[i] - obj.sol(i,0)) > infTop) {
                        infTop = abs(mu.rhs[i] - obj.sol(i,0));
                    }
                    if (abs(mu.rhs[i]) > infBot) {
                        infBot = abs(mu.rhs[i]);
                    }
                }
                f.close();
                cout << "||X_mumps - X_cim||_inf / ||X_mumps||_inf  = " 
                    << infTop/infBot << endl;
                //for(int i =0; i < mu.n; i++){
                    //cout << scientific << mu.rhs[i] << "\t" << obj.sol(i,0) << endl;
                //}
            }
            mu.job = -2;
            dmumps_c(&mu);
        }

    } else {
        int testMumps;
        mpi::broadcast(world, testMumps, 0);
        MUMPS mu;
        if(testMumps > 0) {
            mu.comm_fortran = MPI_Comm_c2f((MPI_Comm) world);
            mu.par = 1;
            mu.job = -1;
            dmumps_c(&mu);

            if (testMumps == 6)
                mu.job = 6;
            else
                mu.job = 4;
            dmumps_c(&mu);
        }
        try {
            obj.bc(-1);
            obj.bc(1);
            obj.bc(2);
            obj.bc(3);
        } catch(int e) {
            cout << world.rank() << " Error code : " << e << endl;
            if(world.rank() == 6) 
            exit(0);
        }

        if(testMumps) {
            mu.job = -2;
            dmumps_c(&mu);
        }
    }
    world.barrier();


    return 0;
}
