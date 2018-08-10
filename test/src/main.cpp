/**
 * \file main.cpp
 * \brief Test program to launch ABCD with a configuration file.
 * \author M. Zenadi
 * \version 1.0
 * \date August 7th 2018
 *
 * This program launch ABCD with the problem and parameters specified in a
 * configuration file (first argument or "config_file.info" by default)
 *
 */
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

/**
 * \fn int main ()
 * \brief Main function
 *
 * \param 1: configuration file (Optional)
 * \param 2: matrix file (Optional)
 * \param 3: number of partitions (Optional)
 * \param 3: block size for the conjugate gradient (Optional)
 * \return EXIT_SUCCESS - Arrêt normal du programme.
 */
int main(int argc, char* argv[])
{
    int error = 0;

    mpi::environment env(argc, argv);
    mpi::communicator world;

    abcd obj;
////////////////////////////////////////////////////////////////////////////////////////
//	ROOT'S JOB
////////////////////////////////////////////////////////////////////////////////////////
    if(world.rank() == 0) {
        ////////////////////////////////////////////////////////////////////////////////
        //	DISPLAY MPI-OPENMP INFO
        ////////////////////////////////////////////////////////////////////////////////
	// Info MPI
	cout << "Number of MPI Procs: " << world.size() << endl;
        // Info OpenMP
	int tid;
	int nthreads;
        #pragma omp parallel private(tid)
        {
                nthreads = omp_get_num_threads();
		tid = omp_get_thread_num();
		#pragma omp master
                cout << "Number of OpenMP threads: " << nthreads << endl;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //	REAGING CONFIG_FILE
        ////////////////////////////////////////////////////////////////////////////////
	/* test existence */
        string config_file;
        if(argc < 2) {
            config_file = "config_file.info";
        } else if(argc >= 2){
            config_file = argv[1];
        }

        /* open and parse the config file */
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

        ////////////////////////////////////////////////////////////////////////////////
        //	GENERAL
        ////////////////////////////////////////////////////////////////////////////////
        obj.write_problem   = pt.get<string>("write_problem", "");
        obj.icntl[Controls::verbose_level] =  pt.get<int>("verbose_level", 0);
        obj.icntl[Controls::mumps_verbose] =  pt.get<int>("mumps_verbose", 0);
        obj.log_output = pt.get<string>("log_filename", "");


        ////////////////////////////////////////////////////////////////////////////////
        //	REAGING THE MATRIX, RHS AND STARTING POINT
        ////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////
        //	open and read the matrix file
        //////////////////////////////////////////
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

        /* allocate the arrays */
        obj.irn = new int[obj.nz];
        obj.jcn = new int[obj.nz];
        obj.val = new double[obj.nz];
        obj.start_index = 1;

        mm_read_mtx_crd_data(f, obj.m, obj.n, obj.nz, obj.irn, obj.jcn, obj.val, mat_code);

        /* Filtering zero-valued elements */
        int newnnz=0;
        for(int i=0;i<obj.nz;i++){
                if(obj.val[i] != 0){
                        obj.irn[newnnz] = obj.irn[i];
                        obj.jcn[newnnz] = obj.jcn[i];
                        obj.val[newnnz] = obj.val[i];
                        newnnz++;
                }
        }

        /* reallocate the arrays ? */
        obj.nz = newnnz;
        obj.irn = (int *) realloc(obj.irn, obj.nz * sizeof(int));
        obj.jcn = (int *) realloc(obj.jcn, obj.nz * sizeof(int));
        obj.val = (double *) realloc(obj.val, obj.nz * sizeof(double));

        cout << "Matrix information : ";
        cout << matrix_file <<" m = " << obj.m << "; n = " << obj.n << "; nz = " << obj.nz << endl;

        fclose(f);

        //////////////////////////////////////////
        //	RHS_FILE
        //////////////////////////////////////////
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
            int ret;
            for(int i = 0; i < m_v*nb_v; i++){
                ret = fscanf(rhs_f, "%lf", &cv);
                obj.rhs[i] = cv;
            }
            fclose(rhs_f);
        }

        //////////////////////////////////////////
        //	START_FILE
        //////////////////////////////////////////
        boost::optional<string> start_file = pt.get_optional<string>("system.start_file");
	if(start_file){
            FILE *strt_f = fopen(start_file->c_str(), "r");

            if(strt_f == NULL){
                cerr << "Error opening the file '"<< *start_file << "'" << endl;
                exit(-1);
            }

            MM_typecode rhs_code;
            int nb_v, m_v;

            mm_read_banner(strt_f, &rhs_code);
            mm_read_mtx_array_size(strt_f, &m_v, &nb_v);

	    obj.Xk = MV_ColMat_double(m_v, nb_v, 0);
            double *xkptr = obj.Xk.ptr();

            cout << "Reading "<< m_v << " values for each of the " << nb_v << " start vector" << endl;

            double cv;
            int ret;
            for(int i = 0; i < m_v*nb_v; i++){
                ret = fscanf(strt_f, "%lf", &cv);
                xkptr[i] = cv;
            }
            fclose(strt_f);
            obj.use_xk=true;
         }

        //////////////////////////////////////////
        //	SOL/BWD_ERR/SCALED_RES
        //////////////////////////////////////////
        boost::optional<string> sol_file = pt.get_optional<string>("system.sol_file");
        boost::optional<string> conv_backward_file = pt.get_optional<string>("system.backward_err_file");
        boost::optional<string> conv_scaled_file = pt.get_optional<string>("system.scaled_residual_file");

        ////////////////////////////////////////////////////////////////////////////////
        //	TEST SOLVING THE PROBLEM WITH MUMPS DIRECTLY
        ////////////////////////////////////////////////////////////////////////////////
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

        ////////////////////////////////////////////////////////////////////////////////
        //	PARTITIONING PART
        ////////////////////////////////////////////////////////////////////////////////
        try{
            ptree::key_type partitioning = pt.get<ptree::key_type>("partitioning");
        } catch (ptree_bad_path e){
            clog << "Error parsing the file, you have to give the partitioning as in the example file" << endl;
            clog << "The what() : " << e.what() << endl << endl;
            exit(-1);
        }

        //////////////////////////////////////////
        //	PARTITIONS
        //////////////////////////////////////////
        obj.icntl[Controls::part_type] = pt.get<int>("partitioning.part_type", 2);

        obj.icntl[Controls::part_guess] = pt.get<int>("partitioning.part_guess", 0);

        if(argc <= 3) obj.icntl[Controls::nbparts] = pt.get<int>("partitioning.nbparts");
        else obj.icntl[Controls::nbparts] = atoi(argv[3]);
        if(obj.icntl[Controls::nbparts] < 0 && obj.icntl[Controls::nbparts] <= obj.m) {
            clog << "Error parsing the file, the number of partitions has to be positive and smaller than the number of rows" << endl;
            exit(-1);
        }
        obj.dcntl[Controls::part_imbalance] = pt.get<double>("partitioning.part_imbalance", 0.5);

        //////////////////////////////////////////
        //	SPECIAL PARTITIONING TYPES
        //////////////////////////////////////////
        /* Manual partitioning */
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
        /* Partitioning file */
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

        //////////////////////////////////////////
        //	FEATURES
        //////////////////////////////////////////
        /* Overlaps */
	obj.icntl[Controls::num_overlap]  = pt.get<int>("partitioning.num_overlap",0);
        if(obj.icntl[Controls::num_overlap] < 0){
            clog << "Error parsing the file, num_overlap must be a positive integer." << endl;
            clog << "Be careful not to input a huge number of overlapping lines, it should not be higher than the smallest partition." << endl;
            exit(-1);
        } else if (obj.icntl[Controls::num_overlap] > 0)
            cout << "Number of overlapping rows: " << obj.icntl[Controls::num_overlap] << endl;

        /* Communication balancing distribution of partitions */
#ifndef NO_METIS
        obj.icntl[Controls::minCommWeight] = pt.get<int>("partitioning.min_comm_weight", 0);
#else
	obj.icntl[Controls::minCommWeight] = 0;
#endif

        /* Enforce Master-Slave scheme: dist_scheme/slave_tol */
        obj.parallel_cg = pt.get<int>("dist_scheme", obj.icntl[Controls::nbparts] < world.size() ? obj.icntl[Controls::nbparts] : world.size());

        obj.icntl[Controls::slave_tol]    = pt.get<int>("partitioning.slave_tol", 0);
        if(obj.icntl[Controls::slave_tol] < 0 || obj.icntl[Controls::slave_tol] > obj.parallel_cg) {
            clog << "Error parsing the file, slave_tol must be a positive integer inferior strictly to the number of Masters." << endl;
            clog << "The number of Masters is the minimum between the number of MPI processes and the number of partitions." << endl;
            exit(-1);
        }
        obj.parallel_cg -= obj.icntl[Controls::slave_tol];

        /* Define masters */
        obj.icntl[Controls::master_def]    = pt.get<int>("partitioning.master_def", 1);
        if(obj.icntl[Controls::master_def] < 0 || obj.icntl[Controls::master_def] > 1){
            clog << "Error parsing the file, master_def must be one of 0 and 1." << endl;
            exit(-1);
        }

        /* Define slaves */
        obj.icntl[Controls::slave_def]    = pt.get<int>("partitioning.slave_def", 2);
        if(obj.icntl[Controls::slave_def] < 0 || obj.icntl[Controls::slave_def] > 2){
            clog << "Error parsing the file, master_def must be one of 0, 1 and 2." << endl;
            exit(-1);
        }
        if (obj.icntl[Controls::master_def] == 0) obj.icntl[Controls::slave_def] = -1;


        ////////////////////////////////////////////////////////////////////////////////
        //	SCALING
        ////////////////////////////////////////////////////////////////////////////////
        obj.icntl[Controls::scaling]    = pt.get<int>("scaling", 2);
        /* manual type */
	if(obj.icntl[Controls::scaling] == -1){
		// Parse number of iterations in option line man_scaling
 		string man_scaling = pt.get<string>("man_scaling", "");
		std::string delimiter = ":";
		size_t pos = 0;
		std::string token;
		int cursor=0;
		while ((pos = man_scaling.find(delimiter)) != std::string::npos && cursor < 4) {
		    token = man_scaling.substr(0, pos);
		    istringstream (token) >> obj.man_scaling[cursor];
		    man_scaling.erase(0, pos + delimiter.length());
		    ++cursor;
		}
                if(obj.man_scaling[0] < 0 || obj.man_scaling[1] < 0 ||
                        obj.man_scaling[2] < 0 || obj.man_scaling[3] < 0) {
                    clog << "Error parsing the file, the number of iterations in scaling must be positive integers." << endl;
                    exit(-1);
                }
        /* predetermined type */
	} else if (obj.icntl[Controls::scaling] == 1) {
		obj.man_scaling[0] = 5; obj.man_scaling[1] = 20;
		obj.man_scaling[2] = 10; obj.man_scaling[3] = 0;
	} else if (obj.icntl[Controls::scaling] == 2) {
		obj.man_scaling[0] = 10; obj.man_scaling[1] = 20;
		obj.man_scaling[2] = 20; obj.man_scaling[3] = 1;
	} else if (obj.icntl[Controls::scaling] == 0){
		obj.man_scaling[0] = 0; obj.man_scaling[1] = 0;
                obj.man_scaling[2] = 0; obj.man_scaling[3] = 0;
	} else {
		clog << "Error: choice of scaling " << obj.icntl[Controls::scaling] << " not recognized" << endl;
                exit(-1);
	}

        ////////////////////////////////////////////////////////////////////////////////
        //	BLOCK CIMMINO
        ////////////////////////////////////////////////////////////////////////////////
        if(argc <= 4) obj.icntl[Controls::block_size] = pt.get<int>("system.block_size", 1);
        else obj.icntl[Controls::block_size] = atoi(argv[4]);
        cout << "Block Size: " << obj.icntl[Controls::block_size] << endl;

        obj.icntl[Controls::itmax] = pt.get<int>("system.itmax", 2000);
        obj.dcntl[Controls::threshold] = pt.get<double>("system.threshold", 1e-12);

        /* scaling factor on the identity of the augmented sub-systems */
        obj.dcntl[Controls::alpha]    = pt.get<double>("system.alpha", 1.0);

        ////////////////////////////////////////////////////////////////////////////////
        //	AUGMENTATION (ABCD)
        ////////////////////////////////////////////////////////////////////////////////
        boost::optional<ptree::key_type> augmentation = pt.get_optional<ptree::key_type>("augmentation");

        if(augmentation){
            obj.icntl[Controls::aug_type]   = pt.get<int>("augmentation.aug_type", 2);
            obj.icntl[Controls::aug_blocking]   = pt.get<int>("augmentation.aug_blocking", 256);
	    cout << "Augmentation type: "<< obj.icntl[Controls::aug_type] << endl ;
	    if(obj.icntl[Controls::aug_type] > 0) cout << "aug_blocking " << obj.icntl[Controls::aug_blocking] << endl;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //	WORK IN PROGRESS (WIP)
        ////////////////////////////////////////////////////////////////////////////////
#ifdef WIP
        /* Use Generalized Gram-Schmidt instead of Generalized QR for the stabilization of CG */
        obj.icntl[Controls::use_gmgs2]  = pt.get<int>("use_gmgs2", 0);
        /* ???? */
        obj.icntl[Controls::exploit_sparcity]  = pt.get<int>("exploit_sparcity", 0);
        if(augmentation){
            /* Only analysis of the augmentation is needed: program stops after augmentation */
            obj.icntl[Controls::aug_analysis]  = pt.get<int>("augmentation.analysis", 0);
            /*
             * Solve with classic BCG on ABCD system if not zero. Ongoing project with
             * preconditioning of S
             */
            obj.icntl[Controls::aug_project]   = pt.get<int>("augmentation.project_only", 0);
            /* deprecated use of dense right hand side when solving for S */
            obj.icntl[Controls::aug_dense]     = pt.get<int>("augmentation.denserhs", 0);
            /*
             * Solve S in an iterative way
             * 0: S solved with MUMPS
             * 1: S solved iteratively, no hard filtering of columns in C, only array of selected columns needed
             * 2: idem but columns hard filtered
             */
            obj.icntl[Controls::aug_iterative] = pt.get<int>("augmentation.iterative", 0);
            /* Threshold for the filtering of the augmentation C */
            obj.dcntl[Controls::aug_filter]    = pt.get<double>("augmentation.filtering", 0.0);
            /*
             * Threshold for the filtering of S which will be solved iteratively (possibly after
             * preconditioning by a matrix M: ongoing project)
             */
            obj.dcntl[Controls::aug_precond]   = pt.get<double>("augmentation.precond", 0.0);
            obj.write_s                        = pt.get<string>("augmentation.write_s", "");
        }
#endif //WIP


        ////////////////////////////////////////////////////////////////////////////////
        //	LAUNCH ABCD SOLVER
        ////////////////////////////////////////////////////////////////////////////////
        try {
            //////////////////////////////////////
            //	ABCD
            //////////////////////////////////////
            obj(-1);

            double t = MPI_Wtime();
            obj(1);
            obj(2);

            obj.nrhs = 1;

            obj(3);

            clog << "Total time: " << MPI_Wtime() - t << endl;

            //////////////////////////////////////
            //	SOL/BWD/SCALED_RES
            //////////////////////////////////////
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
	/* In case of error in ABCD */
        } catch(std::runtime_error e) {
            cout << world.rank() << " Error code : " << e.what() << endl;
            mpi::broadcast(world, error, 0);
            error = 1;
        }

        //////////////////////////////////////////
        //	COMPARE WITH MUMPS SOL ?
        //////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////////////
//	OTHER'S JOB
////////////////////////////////////////////////////////////////////////////////////////
    } else {
        int testMumps;
        mpi::broadcast(world, testMumps, 0);
        MUMPS mu;
        ////////////////////////////////////////////////////////////////////////////////
        //	TEST SOLVING THE PROBLEM WITH MUMPS DIRECTLY
        ////////////////////////////////////////////////////////////////////////////////
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
        ////////////////////////////////////////////////////////////////////////////////
        //	LAUNCH ABCD SOLVER
        ////////////////////////////////////////////////////////////////////////////////
        try {
            obj(-1);
            obj(1);
            obj(2);
            obj(3);
        } catch(std::runtime_error e) {
            cout << world.rank() << " Error code : " << e.what() << endl;
            error=1;
        }
    }
    world.barrier();

    return error;
}
