#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/exceptions.hpp>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>

#include "abcd.h"
#include "mumps.h"

using namespace std;
using namespace boost::property_tree;


extern "C" {
    #include "mmio.h"
    int mm_read_mtx_crd_size(_IO_FILE*, int*, int*, int*);
    int mm_read_mtx_crd_data(_IO_FILE*, int, int, int, int*, int*, double*, char*);
    int mm_read_banner(_IO_FILE*, char (*) [4]);
    int mm_read_mtx_array_size(_IO_FILE*, int*, int*);
}



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
        obj.start_index = 1;

        mm_read_mtx_crd_data(f, obj.m, obj.n, obj.nz, obj.irn, obj.jcn, obj.val, mat_code);
        
        // Filtering zero-valued elements 
        int newnnz=0;
        for(int i=0;i<obj.nz;i++){
                if(obj.val[i] != 0){
                        obj.irn[newnnz] = obj.irn[i];
                        obj.jcn[newnnz] = obj.jcn[i];
                        obj.val[newnnz] = obj.val[i];
                        newnnz++;
                }
        }

        obj.nz = newnnz;
        obj.irn = (int *) realloc(obj.irn, obj.nz * sizeof(int));
        obj.jcn = (int *) realloc(obj.jcn, obj.nz * sizeof(int));
        obj.val = (double *) realloc(obj.val, obj.nz * sizeof(double));

        cout << "Matrix information : ";
        cout << matrix_file <<" m = " << obj.m << "; n = " << obj.n << "; nz = " << obj.nz << endl;

        fclose(f);

        //read the rhs here!
        boost::optional<string> rhs_file = pt.get_optional<string>("system.rhs_file");
        boost::optional<string> sol_file = pt.get_optional<string>("system.sol_file");
        boost::optional<string> conv_backward_file = pt.get_optional<string>("system.backward_err_file");
        boost::optional<string> conv_scaled_file = pt.get_optional<string>("system.scaled_residual_file");
        
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
            int ret;
            for(int i = 0; i < m_v*nb_v; i++){
                ret = fscanf(rhs_f, "%lf", &cv);
                obj.rhs[i] = cv;
            }
            //
            fclose(rhs_f);
        }

        int testMumps =(int) pt.get<bool>("test_mumps", false);
        double minMumps = pt.get<double>("min_mumps", 0);
        double *mumps_rhs;
        int mumps_n = obj.n;
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

            mu.setIcntl(1, 0);
            mu.setIcntl(2, 0);
            mu.setIcntl(3, 0);

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

            clog << "=============================" << endl;
            clog << "LAUNCHING MUMPS ON THE MATRIX" << endl;
            clog << "-----------------------------" << endl;
            double t = MPI_Wtime();
            dmumps_c(&mu);
            clog << "MUMPS RUN FINISHED in " << MPI_Wtime() - t << endl;
            clog << "=============================" << endl;

            if (obj.rhs != NULL){
                mumps_rhs = new double[mu.n];
                std::copy(mu.rhs, mu.rhs + mu.n, mumps_rhs);
            }

            mu.job = -2;
            dmumps_c(&mu);


            if(MPI_Wtime() - t < minMumps && mu.getInfo(1) == 0){
                clog << "MUMPS WAS TOO FAST, RUN!!!!" << endl;
                exit(0);
            }
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

        if(argc <= 3) obj.icntl[Controls::nbparts] = pt.get<int>("partitioning.nbparts");
        else obj.icntl[Controls::nbparts] = atoi(argv[3]);

        if(obj.icntl[Controls::nbparts] < 0 && obj.icntl[Controls::nbparts] <= obj.m) {
            clog << "Error parsing the file, the number of partitions has to be positive and smaller than the number of rows" << endl;
            exit(-1);
        }

        obj.icntl[Controls::part_type] = pt.get<int>("partitioning.part_type", 2);
        obj.icntl[Controls::part_guess] = pt.get<int>("partitioning.part_guess", 0);
        obj.dcntl[Controls::part_imbalance] = pt.get<double>("partitioning.part_imbalance", 0.5);

	if(obj.icntl[Controls::part_type] ==4){
 		string partvector = pt.get<string>("partitioning.partvector", "");;
		obj.partvec   = new int[obj.n];
		const char *cstr = partvector.c_str();
         	std::fstream myfile(cstr, std::ios_base::in);
	        for(int z =0; z < obj.n; z++) {
               		myfile >> obj.partvec[z];
         	}
		myfile.close();		
	}

	cout << "Number of Parts: " << obj.icntl[Controls::nbparts] << " - Partition method: "<< obj.icntl[Controls::part_type] << endl;
	cout << "Number of MPI Procs: " << world.size() << endl;

        if(obj.icntl[Controls::part_type] == 1){
            string parts = pt.get<string>("partitioning.partsfile", "");

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
                    obj.nbrows.push_back(atoi(v.first.data()));
                }

            } else {
                ifstream f;
                f.open(parts.c_str());

                for(unsigned k = 0; k < (unsigned int)obj.icntl[Controls::nbparts]; k++) {
                    string l;
                    getline(f, l);
                    obj.nbrows.push_back(atoi(l.c_str()));
                }

                f.close();
            }

            if(obj.nbrows.size() != (size_t)obj.icntl[Controls::nbparts]){
                clog << "Error parsing the file, nbparts is different from the partitioning description" << endl;
                exit(-1);
            }
        }

        obj.write_problem   = pt.get<string>("write_problem", "");

#ifdef WIP
        obj.icntl[Controls::exploit_sparcity]    = pt.get<int>("exploit_sparcity", 1);
#endif //WIP

        obj.icntl[Controls::scaling]    = pt.get<int>("scaling", 2);

        boost::optional<ptree::key_type> augmentation = pt.get_optional<ptree::key_type>("augmentation");

        if(augmentation){
            obj.icntl[Controls::aug_type]   = pt.get<int>("augmentation.aug_type", 2);
            obj.icntl[Controls::aug_blocking]   = pt.get<int>("augmentation.aug_blocking", 256);
	    cout << "Augmentation type: "<< obj.icntl[Controls::aug_type] << endl ;
	    if(obj.icntl[Controls::aug_type] > 0) cout << "aug_blocking " << obj.icntl[Controls::aug_blocking] << endl;
#ifdef WIP
            obj.icntl[Controls::aug_analysis]   = pt.get<int>("augmentation.analysis", 0);
            obj.dcntl[Controls::aug_filter]   = pt.get<double>("augmentation.filtering", 0.0);
            obj.icntl[Controls::aug_project]   = pt.get<int>("augmentation.project_only", 0);
            obj.icntl[Controls::aug_dense]   = pt.get<int>("augmentation.denserhs", 0);
            obj.icntl[Controls::aug_iterative]   = pt.get<int>("augmentation.iterative", 0);
            obj.dcntl[Controls::aug_precond]   = pt.get<double>("augmentation.precond", 0.0);
            obj.write_s     = pt.get<string>("augmentation.write_s", "");
#endif //WIP
        }

        obj.parallel_cg = pt.get<int>("dist_scheme", obj.icntl[Controls::nbparts] < world.size() ? obj.icntl[Controls::nbparts] : world.size());

        bool error = false;
        try {
            

            obj.icntl[Controls::verbose_level] =  pt.get<int>("verbose_level", 0);
            obj.log_output = pt.get<string>("log_filename", "");

            obj(-1);

            double t = MPI_Wtime();
            obj(1);
            obj(2);
            
            obj.nrhs = 1;

            if(argc <= 4) obj.icntl[Controls::block_size] = pt.get<int>("system.block_size", 1);
            else obj.icntl[Controls::block_size] = atoi(argv[4]);

            obj.icntl[Controls::itmax] = pt.get<int>("system.itmax", 2000);
            obj.dcntl[Controls::threshold] = pt.get<double>("system.threshold", 1e-12);

	    cout << "Block Size: " << obj.icntl[Controls::block_size] << endl;

            // obj.icntl[Controls::verbose] =  pt.get<int>("solve_verbose", 0);
            obj(3);

            clog << "Total time: " << MPI_Wtime() - t << endl;

            if (sol_file) {
                ofstream f; 
                f.open(sol_file->c_str());
                for(int i = 0; i < obj.n_o; i++) {
                    f << obj.sol[i] << "\n";
                }
                f.close();
            }

            if (conv_backward_file) {
                ofstream f;
                f.open(conv_backward_file->c_str());
                f << "Iteration\t BackwardErr\n";
                for(size_t i = 0; i < obj.rhoVector.size() ; i++) {
                    f << i << "\t" << obj.rhoVector[i] << "\n";
                }
                f.close();
            }

            if (conv_scaled_file) {
                ofstream f;
                f.open(conv_scaled_file->c_str());
                f << "Iteration\t Scaled Residual\n";
                for(size_t i = 0; i < obj.scaledResidualVector.size() ; i++) {
                    f << i << "\t" << obj.scaledResidualVector[i] << "\n";
                }
                f.close();
            }

        } catch(std::runtime_error e) {
            cout << world.rank() << " Error code : " << e.what() << endl;
            mpi::broadcast(world, error, 0);
            error = true;
        }

        if(testMumps && !error) {
            double infTop = 0, infBot = 0;
            if(mumps_rhs != NULL){
                ofstream f; 
                f.open("/tmp/out_comp");
                for(int i = 0; i < mumps_n; i++) {
                    f << mumps_rhs[i] << "\t" << obj.sol[i] << "\n";
                    if (abs(mumps_rhs[i] - obj.sol[i]) > infTop) {
                        infTop = abs(mumps_rhs[i] - obj.sol[i]);
                    }
                    if (abs(mumps_rhs[i]) > infBot) {
                        infBot = abs(mumps_rhs[i]);
                    }
                }
                f.close();
                cout << "||X_mumps - X_cim||_inf / ||X_mumps||_inf  = " 
                    << infTop/infBot << endl;
                //for(int i =0; i < mu.n; i++){
                    //cout << scientific << mu.rhs[i] << "\t" << obj.sol(i,0) << endl;
                //}
            }
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
            mu.job = -2;
            dmumps_c(&mu);
        }
        try {
            obj(-1);
            obj(1);
            obj(2);
            obj(3);
        } catch(std::runtime_error e) {
            cout << world.rank() << " Error code : " << e.what() << endl;
        }
    }
    world.barrier();


    return 0;
}
