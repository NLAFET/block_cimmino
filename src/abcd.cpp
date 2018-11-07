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

/*!
 * \file abcd.cpp
 * \brief Main file of the abcd class with gateway function and the main phases
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi, S. Cayrols
 * \version 1.0
 */

#include <abcd.h>
#include <algorithm>
#include <exception>

#include <ctime>
using namespace std;

_INITIALIZE_EASYLOGGINGPP

/*!
 *  \brief Default Constructor of abcd class
 *
 *  Default constructor of the abcd class, sets default values for
 *  attributes.
 *
 */
abcd::abcd()
{
    /* General + write_problem,write_s,log_output */
    last_called_job = -2;
    configure_logger(log_output);

    /* communication */
    mpi::communicator world;
    comm = world;
    parallel_cg = 0; // #masters

    /* Problem */
    start_index = 0;
    // matrix + m,n,nz,sym
    irn = nullptr;
    jcn = nullptr;
    val = nullptr;
    // solution + Xk,sol,rhoVector,scaledResidualVector
    use_xk = false; // CG starting vector specified
    use_xf = false; // ???
    // right hand side
    nrhs = 1;
    rhs = nullptr;
    // augmentation + S_rows,S_cols,S_vals
    size_c = 0;

    /* controls initialization */
    icntl.assign(30, 0);
    dcntl.assign(20, 0);

    // general
    icntl[Controls::verbose_level] = 0;
    icntl[Controls::mumps_verbose] = 0;
    // partitions + strow,nbrows
    icntl[Controls::part_type] = 2;
    icntl[Controls::part_guess] = 1;
    icntl[Controls::nbparts] = 4;
    dcntl[Controls::part_imbalance] = 1.5;
    icntl[Controls::num_overlap] = 0;
    icntl[Controls::overlap_strategy] = 0;
    icntl[Controls::minCommWeight] = 0;
    icntl[Controls::slave_tol] = 0;
    icntl[Controls::master_def] = 0;
    icntl[Controls::slave_def] = -1;
    // scaling
    man_scaling.assign(4, 0);
    icntl[Controls::scaling] = 1;
    man_scaling[0]=5; man_scaling[1]=20;
    man_scaling[2]=10; man_scaling[3]=1;
    // block Cimmino
    icntl[Controls::block_size] = 1;
    icntl[Controls::itmax] = 1000;
    dcntl[Controls::threshold] = 1e-12;
    dcntl[Controls::alpha] = 1.0;
    // augmentation
    icntl[Controls::aug_type] = 0;
    icntl[Controls::aug_blocking] = 256;
#ifdef WIP
    // WIP
    icntl[Controls::exploit_sparcity] = 1;
    icntl[Controls::aug_analysis] = 0;
    icntl[Controls::aug_project] = 0;
    icntl[Controls::aug_dense] = 0
    icntl[Controls::aug_iterative] = 0
    dcntl[Controls::aug_filter] = 0.0;
    dcntl[Controls::aug_precond] = 0.0;
#endif //WIP

    /* infos initialization */
    info.assign(10, 0);
    dinfo.assign(5, 0);
    info[Controls::status] = 0;
    info[Controls::nb_iter] = 0;
}    /* ----- end of constructor abcd::abcd ----- */

/*!
 *  \brief Default Destructor of abcd class
 *
 *  Default destructor of the abcd class, deinitiliaze MUMPS.
 *
 */
abcd::~abcd()
{
  if (mumps.initialized) {
    mumps(-2);
  }
}    /* ----- end of destructor abcd::~abcd ----- */

/*!
 *  \brief Initialize file logger
 *
 *  Initialize the logger to a file if the file is specified by the user
 *  (log_filename)
 *
 */
void abcd::initializeLog() {
        mpi::broadcast(comm, log_output, 0);
        if (comm.rank() > 0 && log_output != "") {
            // find any ".log"
            int ext = log_output.rfind(".log");

            std::string s = boost::lexical_cast<std::string>(comm.rank());

            if (ext != std::string::npos) {
                log_output.replace(ext, 4, "_" + s + ".log");
            } else {
                log_output = log_output + "_" + s;
            }
        }
        logger_set_filename(log_output);
}    /* ----- end of method abcd::initializeLog ----- */

/*!
 *  \brief Creates the internal matrix and rhs from user's data
 *
 *  Centralized on the root:
 *  Initialize the matrix as a Sparselib one based on arrays irn/jcn/val
 *  taking the symmetric/non-symmetric format of the input into account.
 *  If no right hand side was input, also initialize it as A*Xf where
 *  Xf is a m*nrhs matrix with column vectors of values (i+1)/m
 *
 */
int abcd::initializeMatrix()
{
    // Initialization centralized on the root
    if(comm.rank() != 0) return 0;

    LINFO << "*----------------------------------*";
    LINFO << "> Matrix initialization";

    /* Check that the matrix data is present */
    // Allocated matrix ?
    if(irn == nullptr || jcn == nullptr || val == nullptr) {
        // Hey!! where is my data?
        info[Controls::status] = -1;
        LOG_IF(irn == nullptr, ERROR) << "irn is not allocated";
        LOG_IF(jcn == nullptr, ERROR) << "jcn is not allocated";
        LOG_IF(val == nullptr, ERROR) << "val is not allocated";
        throw std::runtime_error("Unallocated matrix vectors");
    }
    LINFO << "M  = " << m << "  N  = " << n << "  NZ = " << nz;

    // Matrix size ok ?
    if(m <= 0 || n <= 0 || nz <= 0){
        info[Controls::status] = -2;
        LOG_IF(m <= 0, ERROR) << "m is negative or zero";
        LOG_IF(n <= 0, ERROR) << "n is negative or zero";
        LOG_IF(nz<= 0, ERROR) << "nz is negative or zero";
        throw std::range_error("Errornous information about the matrix");
    }

    double t = MPI_Wtime();

    LINFO << "Using " << start_index << "-based arrays";

    /// @TODO CHeck that the matrix is not structurally singular
    // only (lower or upper) triangular part provided
    if(sym) {
        // over estimate nz: we should remove diagonal normally but this allows to:
        //  - cope with the zero diagonal elements
        //  - take duplicated entries into account
        int     *t_irn = new int[2*nz];
        int     *t_jcn = new int[2*nz];
        double  *t_val = new double[2*nz];
        int     t_nz = 0;
        for(int k = 0; k < nz; ++k) {
            // transform start_index-indexes into 0-indexes for C++
            irn[k] -= start_index;
            jcn[k] -= start_index;

            t_irn[t_nz] = irn[k];
            t_jcn[t_nz] = jcn[k];
            t_val[t_nz] = val[k];
            t_nz++;
            // duplication of the triangular part to get the full matrix (irn <-> jcn)
            if(irn[k] != jcn[k]){
                t_irn[t_nz] = jcn[k];
                t_jcn[t_nz] = irn[k];
                t_val[t_nz] = val[k];
                t_nz++;
            }
        }
        // Sparselib matrix initialization with actual nz
        nz = t_nz;
        Coord_Mat_double t_A;
        t_A = Coord_Mat_double(m, n, nz, t_val, t_irn, t_jcn, MV_Matrix_::ref);
        A = CompRow_Mat_double(t_A);

        LINFO << "> Matrix initialized in " << setprecision(2) << MPI_Wtime() - t << "s.";
        delete[] t_irn;
        delete[] t_jcn;
        delete[] t_val;
    } else {
        t = MPI_Wtime();
        for(int i=0; i<nz; ++i){
            // transform start_index-indexes into 0-indexes for C++
            irn[i] -= start_index;
            jcn[i] -= start_index;
        }
        Coord_Mat_double t_A;
        t_A = Coord_Mat_double(m, n, nz, val, irn, jcn, MV_Matrix_::ref);
        A = CompRow_Mat_double(t_A);

        LINFO << "> Matrix initialized in " << setprecision(2) << MPI_Wtime() - t << "s.";
    }

    // Save original size of matrix
    n_o = n;
    m_o = m;
    nz_o = nz;

    // If no RHS specified, create an artificial with values (i+1)/m
    if(rhs == nullptr){
        LINFO << "-- No RHS specified, the new one will be created --";
        t = MPI_Wtime();
	// Create a matrix Xf where columns are vectors of values (i+1)/m
        Xf = MV_ColMat_double(m, nrhs);
        for(int j = 0; j < nrhs; j++){
            VECTOR_double xf_col(m);
            for(int i = 0; i < m; i++) {
                xf_col[i] = (double)(i+1) / m;
            }
            Xf.setCol(xf_col, j);
        }
	// Create the RHS as rhs=A*Xf
        MV_ColMat_double BB = smv(A, Xf);
        // To Do : make it for Multiple RHS
        rhs = new double[m * nrhs];
        double *BBref = BB.ptr();
        for(int i = 0; i < m; i++){
            rhs[i] = BBref[i];
        }
/*      FOR MULTIPLE RHS ??
        for(int j=0; j<nrhs; ++j) {
            for(int i = 0; i < m; i++){
                rhs[j*m+i] = BBref[j*m+i];
            }
        }*/
        LINFO << "> RHS initialized in " << setprecision(2) << MPI_Wtime() - t << "s.";
    }
    return 0;
}    /* ----- end of method abcd::initializeMatrix ----- */

/*!
 *  \brief Scales and partitions the matrix then analyses the structure of the partitions
 *
 *  Centralized on the root:
 *  Method to scale the matrix, find its partitioning, create the partitions and augmentation.
 *
 */
int abcd::preprocessMatrix()
{
    LINFO << "*----------------------------------*";
    LINFO << "> Starting Preprocessing            ";
    double t = MPI_Wtime();
    double tot = t;
    int err = 0;

    // Preprocessing centralized on the root
    // + Error checking
    if(comm.rank() != 0) {
        mpi::broadcast(comm, err, 0);
        if (err != 0) {
            info[Controls::status] = err;
            throw std::runtime_error("The master asked me to stop.");
        }
        return err;
    }

    // If #masters not set, min(#parts, #MPI) by default
    if (parallel_cg == 0) {
        parallel_cg = icntl[Controls::nbparts] < comm.size() ? icntl[Controls::nbparts] : comm.size();
    }

    // Scaling
    t = MPI_Wtime();

    abcd::scaling();

    LINFO << "> Time to scale the matrix: "
          << MPI_Wtime() - t << "s.";

    // Find partitioning
    t = MPI_Wtime();

    abcd::partitionMatrix();

    LINFO << "> Time to partition the matrix: "
          << MPI_Wtime() - t << "s.";

    // Create partitions/augmentation
    abcd::analyseFrame();
    LINFO << "> Total time to preprocess: " << MPI_Wtime() - tot << "s.";
    LINFO << "*----------------------------------*";

    // everything is alright, tell the others that we're done here
    mpi::broadcast(comm, err, 0);

    return 0;
}    /* ----- end of method abcd::preprocessMatrix ----- */

/*!
 *  \brief Creates augmented systems and factorizes them
 *
 *  Creates the augmented systems based on the partitions and factorizes them
 *  with MUMPS.
 *
 */
int abcd::factorizeAugmentedSystems()
{
    // Create the group of CG instances
    abcd::createInterCommunicators();

    double t = MPI_Wtime();

    if(instance_type == 0) {
        // Broadcast original size of the matrix
        mpi::broadcast(inter_comm, m_o, 0);
        mpi::broadcast(inter_comm, n_o, 0);
        mpi::broadcast(inter_comm, nz_o, 0);
        // Broadcast controls
        mpi::broadcast(inter_comm, icntl, 0);
        mpi::broadcast(inter_comm, dcntl, 0);
        // Broadcast size of the augmentation
        if( icntl[Controls::aug_type] != 0 )
            mpi::broadcast(inter_comm, size_c, 0);
    }

    /*
     * distribute partitions and slaves to masters then compute interconnections and
     * initialize MUMPS
     */
    if(instance_type == 0) {
        abcd::distributeData();
        abcd::createInterconnections();
    }

    abcd::initializeDirectSolver(icntl[Controls::inner_solver]);

    // Do the Analysis in all cases except when the process is a Master with no slaves and
    // there are slaves
    // SCAYROLS_ADD add test on the solver since the analyse should have already done during
    //              the init of the direct solver
    if (! (instance_type == 0 && comm.size() > parallel_cg && my_slaves.size() == 0) 
        && icntl[Controls::inner_solver] == MUMPS_SOLVER_TYPE) {
        if(inter_comm.rank() == 0 && instance_type == 0) {
            LINFO << "Launching MUMPS analysis";
        }
        abcd::analyseAugmentedSystems(mumps);
    }
    if(IRANK == 0){
        LINFO << "Initialization time : " << MPI_Wtime() - t;
    }

    switch(icntl[Controls::inner_solver]){
      case MUMPS_SOLVER_TYPE :
        if(inter_comm.rank() == 0 && instance_type == 0){
          LINFO << "Launching MUMPS factorization";
        }

        t = MPI_Wtime();
        abcd::factorizeAugmentedSystems(mumps);
        break;
      case SPLDLT_SOLVER_TYPE :
        if(inter_comm.rank() == 0 && instance_type == 0){
          LINFO << "Launching SpLDLT factorization";
        }

        t = MPI_Wtime();
        abcd::factorizeAugmentedSystems(inner_solver);
        break;
    }

    if(IRANK == 0){
      LINFO << "Factorization time : " << MPI_Wtime() - t;
    }

    return 0;
}    /* ----- end of method abcd::factorizeAugmentedSystems ----- */

/*!
 *  \brief Runs either BCG or ABCD solve depending on what we want
 *
 *  Distribute the right hand side then solve the system using BCG or ABCD depending
 *  on the choice (icntl[Controls::aug_type]).
 *
 */
int abcd::solveSystem()
{
    // Check RHS and Block Size positive
    if (nrhs < 1) {
        info[Controls::status] = -9;
        throw std::runtime_error("Number of right-hand sides should be at least one (1)");
    }
    if(icntl[Controls::block_size] < 1) {
        info[Controls::status] = -10;
        throw std::runtime_error("Block size should be at least one (1)");
    }

    // Hard synchronisation between masters
    if(instance_type == 0) inter_comm.barrier();
    if(inter_comm.rank() == 0 && instance_type == 0){
        LINFO << "Launching Solve";
    }

    // Broadcast MUMPS controls to all masters
    if(instance_type == 0) {
        mpi::broadcast(inter_comm, icntl, 0);
        mpi::broadcast(inter_comm, dcntl, 0);

        // If augmentation empty, get back to classic BCG
        if(size_c == 0 && inter_comm.rank() == 0 && icntl[Controls::aug_type] != 0){
            LINFO << "Size of S is 0, therefore launching bcg";
            icntl[Controls::aug_type] = 0;
        }

        // If we are running ABCD, we have to ensure that we have no block-size > 1
        // Block Cimmino
        if(icntl[Controls::aug_type] == 0
#ifdef WIP
           || icntl[Controls::aug_project] != 0
#endif //WIP
                ){
            abcd::distributeRhs();
            abcd::bcg(B);
        // Augmented Block Cimmino
        } else{
            int temp_bs = icntl[Controls::block_size];
            icntl[Controls::block_size] = 1;
            abcd::distributeRhs();
            abcd::solveABCD(B);
            icntl[Controls::block_size] = temp_bs;
        }

        int job = -1;
        mpi::broadcast(intra_comm, job, 0);

        // Display information on the accuracy of the solution
        if(inter_comm.rank() == 0){
            LINFO1 << "Backward error       : " << scientific << dinfo[Controls::backward];
            LINFO1 << "||r||_inf            : " << scientific << dinfo[Controls::residual];
            LINFO1 << "||r||_inf/||b||_inf  : " << scientific << dinfo[Controls::scaled_residual];
            if (Xf.dim(0) != 0)
                LINFO << "Forward error        : " <<
                    scientific << dinfo[Controls::forward_error];
        }
    } else {
        if(icntl[Controls::aug_type] != 0 ){
            int temp_bs = icntl[Controls::block_size];
            icntl[Controls::block_size] = 1;
            abcd::waitForSolve();
            icntl[Controls::block_size] = temp_bs;
        } else {
            abcd::waitForSolve();
        }
    }

    return 0;
}    /* ----- end of method abcd::solveSystem ----- */

/*!
 *  \brief Gateway function
 *
 *  The gateway function that launches all other options.
 *
 *  \param job_id: the parameter choosing which phase of the solver to launch
 *  \return 0 no error occurred
 */
int abcd::operator()(int job_id)
{
    info[Controls::status] = 0;

    LDEBUG3 << "MPI-Process " << comm.rank() << " called job_id = " << job_id;

    ////////////////////////////////////////////////////////////////////////////////////
    //      CHECK RIGHT ORDER OF JOB PARAMETER CALLED
    ////////////////////////////////////////////////////////////////////////////////////
    if ( job_id == 1 && last_called_job != -1 ) {
        info[Controls::status] = -2;
        throw std::runtime_error("Did you forget to call job_id = -1? ");
    }
    if ( job_id == 2 && job_id <= 3 && last_called_job != 1) {
        info[Controls::status] = -2;
        throw std::runtime_error("Did you forget to call job_id = 1? ");
    }
    if ( job_id == 3 && last_called_job != 2 && last_called_job != 3 ) {
        info[Controls::status] = -2;
        throw std::runtime_error("Did you forget to call job_id = 2? ");
    }
    if ( job_id < last_called_job && job_id != 3) {
        info[Controls::status] = -2;
        throw std::runtime_error("You cannot go back in time");
    }
    if ( job_id == last_called_job && job_id != 3) {
        info[Controls::status] = -2;
        throw std::runtime_error("You cannot re-run the same job unless you want to re-run job_id = 3 ");
    }

    ////////////////////////////////////////////////////////////////////////////////////
    //      CALL THE FUNCTION CORRESPONDING TO JOB
    ////////////////////////////////////////////////////////////////////////////////////
    switch(job_id) {
    /* Initialization */
    case -1:
        // in case the user decided to change the log_file name or
        // decided to disable the logging
        initializeLog();

        initializeMatrix();

        // if everything went alright, remember the job
        last_called_job = job_id;

        break;
    /* Preprocessing */
    case 1:
        preprocessMatrix();

        // if everything went alright, remember the job
        last_called_job = job_id;

        break;
    /* Factorization */
    case 2:
        factorizeAugmentedSystems();

        // if everything went alright, remember the job
        last_called_job = job_id;

        break;
    /* Solve */
    case 3:
        solveSystem();

        // if everything went alright, remember the job
        last_called_job = job_id;

        break;
    /* Preprocessing + Factorization */
    case 4:
        // this ensures that we get consistent results
        operator()(1);
        operator()(2);

        // if everything went alright, remember the job
        last_called_job = 2;

        break;
    /* Factorization + Solve */
    case 5:
        operator()(2);
        operator()(3);

        // if everything went alright, remember the job
        last_called_job = 3;

        break;
    /* Preprocessing + Factorization + Solve */
    case 6:
        operator()(1);
        operator()(2);
        operator()(3);

        // if everything went alright, remember the job
        last_called_job = 3;

        break;
    /* Unrecognized job parameter */
    default:
        // Wrong job id
        info[Controls::status] = -1;
        throw std::runtime_error("Wrong job_id id.");
    }
    return info[Controls::status];
}    /* ----- end of method abcd::operator() ----- */
