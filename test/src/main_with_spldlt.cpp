/**
 * \file main.cpp
 * \brief Test program to launch ABCD with a configuration file.
 * \author S. Cayrols
 * \version 1.0
 * \date October 18th 2018
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

  abcd solver;
  //////////////////////////////////////////////////////////////////////////////
  //	ROOT'S JOB
  //////////////////////////////////////////////////////////////////////////////
  if(world.rank() == 0) {
    // Info MPI
    cout << "Number of MPI Procs: " << world.size() << endl;
    // Info OpenMP
    int nthreads;
    #pragma omp single
    {
      nthreads = omp_get_num_threads();
      cout << "Number of OpenMP threads: " << nthreads << endl;
    }

    ////////////////////////////////////////////////////////////////////////////
    //	READING CONFIG_FILE
    ////////////////////////////////////////////////////////////////////////////
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

    ////////////////////////////////////////////////////////////////////////////
    //	GENERAL
    ////////////////////////////////////////////////////////////////////////////
    solver.write_problem                  = pt.get<string>("write_problem", "");
    solver.icntl[Controls::verbose_level] = pt.get<int>("verbose_level", 0);
    solver.icntl[Controls::mumps_verbose] = pt.get<int>("mumps_verbose", 0);
    solver.log_output                     = pt.get<string>("log_filename", "");


    ////////////////////////////////////////////////////////////////////////////
    //	READING THE MATRIX, RHS AND STARTING POINT
    ////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////
    //	open and read the matrix file
    //////////////////////////////////////////
    string matrix_file;
    try {
      if(argc <= 2){
        matrix_file = pt.get<string>("system.matrix_file");
      }
      else{
        matrix_file = argv[2];
      }
    } catch (ptree_bad_path e) {
      clog << "Error parsing the file, you have to give"
        "the matrix name as in the example file" << endl;
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
    mm_read_mtx_crd_size(f, (int *)&solver.m,
        (int *)&solver.n, (int *)&solver.nz);

    /* allocate the arrays */
    solver.irn = new int[solver.nz];
    solver.jcn = new int[solver.nz];
    solver.val = new double[solver.nz];
    solver.start_index = 1;

    mm_read_mtx_crd_data(f, solver.m, solver.n, solver.nz,
        solver.irn, solver.jcn, solver.val, mat_code);

    /* Filtering zero-valued elements */
    int newnnz=0;
    for(int i=0;i<solver.nz;i++){
      if(solver.val[i] != 0){
        solver.irn[newnnz] = solver.irn[i];
        solver.jcn[newnnz] = solver.jcn[i];
        solver.val[newnnz] = solver.val[i];
        newnnz++;
      }
    }

    /* reallocate the arrays ? */
    solver.nz = newnnz;
    solver.irn = (int *) realloc(solver.irn, solver.nz * sizeof(int));
    solver.jcn = (int *) realloc(solver.jcn, solver.nz * sizeof(int));
    solver.val = (double *) realloc(solver.val, solver.nz * sizeof(double));

    cout << "Matrix information : ";
    cout << matrix_file <<" m = " << solver.m << "; n = " << solver.n <<
      "; nz = " << solver.nz << endl;

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

      solver.rhs = new double[nb_v * m_v];
      cout << "Reading "<< m_v << " values for each of the " <<
        nb_v << " rhs" << endl;

      double cv;
      int ret;
      for(int i = 0; i < m_v*nb_v; i++){
        ret = fscanf(rhs_f, "%lf", &cv);
        solver.rhs[i] = cv;
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

      solver.Xk = MV_ColMat_double(m_v, nb_v, 0);
      double *xkptr = solver.Xk.ptr();

      cout << "Reading "<< m_v << " values for each of the " << nb_v << " start vector" << endl;

      double cv;
      int ret;
      for(int i = 0; i < m_v*nb_v; i++){
        ret = fscanf(strt_f, "%lf", &cv);
        xkptr[i] = cv;
      }
      fclose(strt_f);
      solver.use_xk=true;
    }

    //////////////////////////////////////////
    //	SOL/BWD_ERR/SCALED_RES
    //////////////////////////////////////////
    boost::optional<string> sol_file = pt.get_optional<string>("system.sol_file");
    boost::optional<string> conv_backward_file = pt.get_optional<string>("system.backward_err_file");
    boost::optional<string> conv_scaled_file = pt.get_optional<string>("system.scaled_residual_file");

    ////////////////////////////////////////////////////////////////////////////////
    //	PARTITIONING PART
    ////////////////////////////////////////////////////////////////////////////////
    try{
      ptree::key_type partitioning = pt.get<ptree::key_type>("partitioning");
    } catch (ptree_bad_path e){
      clog << "Error parsing the file, you have to give the partitioning"
        "as in the example file" << endl;
      clog << "The what() : " << e.what() << endl << endl;
      exit(-1);
    }

    //////////////////////////////////////////
    //	PARTITIONS
    //////////////////////////////////////////
    solver.icntl[Controls::part_type] = pt.get<int>("partitioning.part_type", 2);

    solver.icntl[Controls::part_guess] = pt.get<int>("partitioning.part_guess", 0);

    if(argc <= 3) solver.icntl[Controls::nbparts] = pt.get<int>("partitioning.nbparts");
    else solver.icntl[Controls::nbparts] = atoi(argv[3]);
    if(solver.icntl[Controls::nbparts] < 0 && solver.icntl[Controls::nbparts] <= solver.m) {
      clog << "Error parsing the file, the number of partitions has to be"
        "positive and smaller than the number of rows" << endl;
      exit(-1);
    }
    solver.dcntl[Controls::part_imbalance] = pt.get<double>("partitioning.part_imbalance", 0.5);

    //////////////////////////////////////////
    //	SPECIAL PARTITIONING TYPES
    //////////////////////////////////////////
    /* Manual partitioning */
    if(solver.icntl[Controls::part_type] == 1){
      string parts = pt.get<string>("partitioning.partsfile", "");

      if(parts.length() == 0){
        try{
          ptree::key_type partitioning = pt.get<ptree::key_type>("partitioning.nbrows");
        } catch (ptree_bad_path e){
          clog << "Error parsing the file, you have to give the number of rows"
            "per partition as in the example file" << endl;
          clog << "The what() : " << e.what() << endl << endl;
          exit(-1);
        }

        BOOST_FOREACH( ptree::value_type v, pt.get_child("partitioning.nbrows") )
        {
          solver.nbrows.push_back(atoi(v.first.data()));
        }

      } else {
        ifstream f;
        f.open(parts.c_str());

        for(unsigned k = 0; k < (unsigned int)solver.icntl[Controls::nbparts]; k++) {
          string l;
          getline(f, l);
          solver.nbrows.push_back(atoi(l.c_str()));
        }

        f.close();
      }

      if(solver.nbrows.size() != (size_t)solver.icntl[Controls::nbparts]){
        clog << "Error parsing the file, nbparts is different from the"
          "partitioning description" << endl;
        exit(-1);
      }
    }
    /* Partitioning file */
    if(solver.icntl[Controls::part_type] ==4){
      string partvector = pt.get<string>("partitioning.partvector", "");;
      solver.partvec   = new int[solver.n];
      const char *cstr = partvector.c_str();
      std::fstream myfile(cstr, std::ios_base::in);
      for(int z =0; z < solver.n; z++) {
        myfile >> solver.partvec[z];
      }
      myfile.close();
    }

    cout << "Number of Parts: " << solver.icntl[Controls::nbparts] <<
      " - Partition method: "<< solver.icntl[Controls::part_type] << endl;

    //////////////////////////////////////////
    //	FEATURES
    //////////////////////////////////////////
    /* Overlaps */
    solver.icntl[Controls::num_overlap]  = pt.get<int>("partitioning.num_overlap",0);
    if(solver.icntl[Controls::num_overlap] < 0){
      clog << "Error parsing the file, num_overlap must be"
        "a positive integer." << endl;
      clog << "Be careful not to input a huge number of overlapping lines,"
           "it should not be higher than the smallest partition." << endl;
      exit(-1);
    } else if(solver.icntl[Controls::num_overlap] > solver.m){
      clog << "Error parsing the file, num_overlap must be less than"
        "the matrix size." << endl;
      clog << "Be careful not to input a huge number of overlapping lines,"
           "it should not be higher than the smallest partition." << endl;
      exit(-1);
    } else if (solver.icntl[Controls::num_overlap] > 0)
      cout << "Number of overlapping rows: " <<
        solver.icntl[Controls::num_overlap] << endl;

    /* Communication balancing distribution of partitions */
#ifndef NO_METIS
    solver.icntl[Controls::minCommWeight] = pt.get<int>("partitioning.min_comm_weight", 0);
#else
	  solver.icntl[Controls::minCommWeight] = 0;
#endif

    /* Enforce Master-Slave scheme: dist_scheme/slave_tol */
    solver.parallel_cg = pt.get<int>("dist_scheme", 
        solver.icntl[Controls::nbparts] < world.size() ?
        solver.icntl[Controls::nbparts] : world.size());

    solver.icntl[Controls::slave_tol]    = pt.get<int>("partitioning.slave_tol", 0);
    if(solver.icntl[Controls::slave_tol] < 0 ||
        solver.icntl[Controls::slave_tol] > solver.parallel_cg) {
      clog << "Error parsing the file, slave_tol must be a positive integer"
        "inferior strictly to the number of Masters." << endl;
      clog << "The number of Masters is the minimum between the number of"
        "MPI processes and the number of partitions." << endl;
      exit(-1);
    }
    solver.parallel_cg -= solver.icntl[Controls::slave_tol];

    /* Define masters */
    solver.icntl[Controls::master_def]    = pt.get<int>("partitioning.master_def", 1);
    if(solver.icntl[Controls::master_def] < 0 || solver.icntl[Controls::master_def] > 1){
      clog << "Error parsing the file, master_def must be"
        "one of 0 and 1." << endl;
      exit(-1);
    }

    /* Define slaves */
    solver.icntl[Controls::slave_def]    = pt.get<int>("partitioning.slave_def", 2);
    if(solver.icntl[Controls::slave_def] < 0 || solver.icntl[Controls::slave_def] > 2){
      clog << "Error parsing the file, master_def must be"
        "one of 0, 1 and 2." << endl;
      exit(-1);
    }
    if (solver.icntl[Controls::master_def] == 0) solver.icntl[Controls::slave_def] = -1;


    ////////////////////////////////////////////////////////////////////////////
    //	SCALING
    ////////////////////////////////////////////////////////////////////////////
    solver.icntl[Controls::scaling]    = pt.get<int>("scaling", 2);
    /* manual type */
    if(solver.icntl[Controls::scaling] == -1){
      // Parse number of iterations in option line man_scaling
      string man_scaling = pt.get<string>("man_scaling", "");
      std::string delimiter = ":";
      size_t pos = 0;
      std::string token;
      int cursor=0;
      while ((pos = man_scaling.find(delimiter)) != std::string::npos && cursor < 4) {
        token = man_scaling.substr(0, pos);
        istringstream (token) >> solver.man_scaling[cursor];
        man_scaling.erase(0, pos + delimiter.length());
        ++cursor;
      }
      if(solver.man_scaling[0] < 0 || solver.man_scaling[1] < 0 ||
          solver.man_scaling[2] < 0 || solver.man_scaling[3] < 0) {
        clog << "Error parsing the file, the number of iterations in"
          "scaling must be positive integers." << endl;
        exit(-1);
      }
      /* predetermined type */
    } else if (solver.icntl[Controls::scaling] == 1) {
      solver.man_scaling[0] = 5; solver.man_scaling[1] = 20;
      solver.man_scaling[2] = 10; solver.man_scaling[3] = 0;
    } else if (solver.icntl[Controls::scaling] == 2) {
      solver.man_scaling[0] = 10; solver.man_scaling[1] = 20;
      solver.man_scaling[2] = 20; solver.man_scaling[3] = 1;
    } else if (solver.icntl[Controls::scaling] == 0){
      solver.man_scaling[0] = 0; solver.man_scaling[1] = 0;
      solver.man_scaling[2] = 0; solver.man_scaling[3] = 0;
    } else {
      clog << "Error: choice of scaling " <<
        solver.icntl[Controls::scaling] << " not recognized" << endl;
      exit(-1);
    }

    ////////////////////////////////////////////////////////////////////////////
    //	BLOCK CIMMINO
    ////////////////////////////////////////////////////////////////////////////
    if(start_file){
      cout << "Starting point for CG specified,"
           "Block Size is changed to 1.\n";
      solver.icntl[Controls::block_size] = 1;
    } else {
      if(argc <= 4) {
        solver.icntl[Controls::block_size] = pt.get<int>("system.block_size", 1);
      } else {
        solver.icntl[Controls::block_size] = atoi(argv[4]);
      }
      cout << "Block Size: " << solver.icntl[Controls::block_size] << endl;
    }

    solver.icntl[Controls::itmax] = pt.get<int>("system.itmax", 2000);
    solver.dcntl[Controls::threshold] = pt.get<double>("system.threshold", 1e-12);

    if(pt.get<int>("system.innerSolver", MUMPS_SOLVER_TYPE) == SPLDLT_SOLVER_TYPE)
      solver.icntl[Controls::innerSolver] = SPLDLT_SOLVER_TYPE;
    else
      solver.icntl[Controls::innerSolver] = MUMPS_SOLVER_TYPE;

    /* scaling factor on the identity of the augmented sub-systems */
    solver.dcntl[Controls::alpha] = pt.get<double>("system.alpha", 1.0);

    ////////////////////////////////////////////////////////////////////////////////
    //	AUGMENTATION (ABCD)
    ////////////////////////////////////////////////////////////////////////////////
    boost::optional<ptree::key_type> augmentation = pt.get_optional<ptree::key_type>("augmentation");

    if(augmentation){
      cout << "ABCD version "<< endl;
      solver.icntl[Controls::aug_type]     = pt.get<int>("augmentation.aug_type", 2);
      solver.icntl[Controls::aug_blocking] = pt.get<int>("augmentation.aug_blocking", 256);
      cout << "Augmentation type: "<< solver.icntl[Controls::aug_type] << endl ;
      if(solver.icntl[Controls::aug_type] > 0) {
        cout << "aug_blocking " << solver.icntl[Controls::aug_blocking] << endl;
      }
    }
    else{
      cout << "Block Cimmino simple"<< endl;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //	LAUNCH ABCD SOLVER
    ////////////////////////////////////////////////////////////////////////////////
    try {
      //////////////////////////////////////
      //	ABCD
      //////////////////////////////////////
      solver(-1);

      double t = MPI_Wtime();
      solver(1);
      solver(2);

      solver.nrhs = 1;

      solver(3);

      clog << "Total time: " << MPI_Wtime() - t << endl;

      //////////////////////////////////////
      //	SOL/BWD/SCALED_RES
      //////////////////////////////////////
      if (sol_file) {
        ofstream f;
        f.open(sol_file->c_str());
        for(int i = 0; i < solver.n_o; i++) {
          f << solver.sol[i] << "\n";
        }
        f.close();
      }

      if (conv_backward_file) {
        ofstream f;
        f.open(conv_backward_file->c_str());
        f << "Iteration\t BackwardErr\n";
        for(size_t i = 0; i < solver.rhoVector.size() ; i++) {
          f << i << "\t" << solver.rhoVector[i] << "\n";
        }
        f.close();
      }

      if (conv_scaled_file) {
        ofstream f;
        f.open(conv_scaled_file->c_str());
        f << "Iteration\t Scaled Residual\n";
        for(size_t i = 0; i < solver.scaledResidualVector.size() ; i++) {
          f << i << "\t" << solver.scaledResidualVector[i] << "\n";
        }
        f.close();
      }
      /* In case of error in ABCD */
    } catch(std::runtime_error e) {
      cout << world.rank() << " Error code : " << e.what() << endl;
      world.abort(1);
      mpi::broadcast(world, error, 0);
      error = 1;
    }

  }
////////////////////////////////////////////////////////////////////////////////
//	OTHER'S JOB
////////////////////////////////////////////////////////////////////////////////
  else {
    ////////////////////////////////////////////////////////////////////////////
    //	LAUNCH ABCD SOLVER
    ////////////////////////////////////////////////////////////////////////////
    try {
      solver(-1);
      solver(1);
      solver(2);
      solver(3);
    } catch(std::runtime_error e) {
      cout << world.rank() << " Error code : " << e.what() << endl;
      world.abort(1);
      mpi::broadcast(world, error, 0);
      error=1;
    }
  }
  world.barrier();

  return error;
}
