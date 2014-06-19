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
#include <algorithm>
#include <exception>

#include <ctime>
using namespace std;

_INITIALIZE_EASYLOGGINGPP

/// Set the defaults in the constructor
abcd::abcd()
{
    last_called_job = -2;
    
    nrhs = 1;
    start_index = 0;
    use_xk = false;
    use_xf = false;
    rhs = nullptr;
    size_c = 0;
    verbose = false;
    runSolveS = false;
    parallel_cg = 0;

    irn = nullptr;
    jcn = nullptr;
    val = nullptr;

    icntl.assign(20, 0);
    dcntl.assign(20, 0);
    info.assign(10, 0);
    dinfo.assign(5, 0);

    icntl[Controls::aug_blocking] = 256;

    icntl[Controls::nbparts] = 4;
    icntl[Controls::part_type] = 2;
    dcntl[Controls::part_imbalance] = 1.5;
    icntl[Controls::part_guess] = 1;
    icntl[Controls::scaling] = 2;
    icntl[Controls::itmax] = 1000;
    icntl[Controls::block_size] = 1;
    dcntl[Controls::threshold] = 1e-12;

    icntl[Controls::verbose_level] = 0;
    info[Controls::status] = 0;
    info[Controls::nb_iter] = 0;

    mpi::communicator world;
    comm = world;


    // Prepare for Logging!
    if (comm.rank() == 0 && log_output == "") {
        time_t t = time(0);
        struct tm * now = localtime(&t);

        char timeString[15];  // space for "HH_MM_SS_dd_mm\0"

        strftime(timeString, sizeof(timeString), "%H_%M_%S_%d_%m", now);

        std::string tt(timeString);
        log_output = "/tmp/log_abcd_" + tt + ".log";
    }
    configure_logger(log_output);
}

abcd::~abcd()
{
  if (mumps.initialized) {
    mumps(-2);
  }
}

/// Creates the internal matrix from user's data
int abcd::initializeMatrix()
{
    LINFO << "*----------------------------------*";
    LINFO << "> Local matrix initialization";
    
    if(comm.rank() != 0) return 0;
    
    // Check that the matrix data is present
    if(irn == nullptr || jcn == nullptr || val == nullptr) {
        // Hey!! where is my data?
        info[Controls::status] = -1;
        LOG_IF(irn == nullptr, ERROR) << "irn is not allocated";
        LOG_IF(jcn == nullptr, ERROR) << "jcn is not allocated";
        LOG_IF(val == nullptr, ERROR) << "val is not allocated";
        throw std::runtime_error("Unallocated matrix vectors");
    }
    LINFO << "M  = " << m << "  N  = " << n << "  NZ = " << nz;

    if(m <= 0 || n <= 0 || nz <= 0){
        info[Controls::status] = -2;
        LOG_IF(m <= 0, ERROR) << "m is negative or zero";
        LOG_IF(n <= 0, ERROR) << "n is negative or zero";
        LOG_IF(nz<= 0, ERROR) << "nz is negative or zero";
        throw std::range_error("Errornous information about the matrix");
    }

    double t = MPI_Wtime();
    int st_point = 1;
    // chose 0-based or 1-based
    if ( irn[0] == 0 || jcn[0] == 0 ) st_point = 0;

    LINFO << "Using " << st_point << "-based arrays";

    /// @TODO CHeck that the matrix is not structurally singular
    if(sym) {
        //over estimate nz
        int     *t_irn = new int[2*nz];
        int     *t_jcn = new int[2*nz];
        double  *t_val = new double[2*nz];
        int     t_nz = 0;
        for(int k = 0; k < nz; ++k) {
            irn[k] -= st_point;
            jcn[k] -= st_point;

            t_irn[t_nz] = irn[k];
            t_jcn[t_nz] = jcn[k];
            t_val[t_nz] = val[k];
            t_nz++;
            if(irn[k] != jcn[k]){
                t_irn[t_nz] = jcn[k];
                t_jcn[t_nz] = irn[k];
                t_val[t_nz] = val[k];
                t_nz++;
            }
        }
        nz = t_nz;
        Coord_Mat_double t_A;
        t_A = Coord_Mat_double(m, n, t_nz, t_val, t_irn, t_jcn);
        t = MPI_Wtime();
        A = CompRow_Mat_double(t_A);

        LINFO << "> Local matrix initialized in " << setprecision(2) << MPI_Wtime() - t << "s.";
        delete[] t_irn;
        delete[] t_jcn;
        delete[] t_val;
    } else {
        t = MPI_Wtime();
        for(int i=0; i<nz; ++i){
            irn[i] -= st_point;
            jcn[i] -= st_point;

        }
        Coord_Mat_double t_A;
        t_A = Coord_Mat_double(m, n, nz, val, irn, jcn, MV_Matrix_::ref);
        A = CompRow_Mat_double(t_A);

        LINFO << "> Local matrix initialized in " << setprecision(2) << MPI_Wtime() - t << "s.";
    }
    
    n_o = n;
    m_o = m;
    nz_o = nz;
        
    return 0; 
}

/// Scales, partitions and analyses the structure of partitions
int abcd::preprocessMatrix()
{

    LINFO << "*----------------------------------*";
    LINFO << "> Starting Preprocessing            ";
    double t = MPI_Wtime();
    double tot = t;
    int err = 0;

    if(comm.rank() != 0) {
        mpi::broadcast(comm, err, 0);
        if (err != 0) {
            info[Controls::status] = err;
            throw std::runtime_error("The master asked me to stop.");
        }
        return err;
    }

    if (parallel_cg == 0) {
        parallel_cg = icntl[Controls::nbparts] < comm.size() ? icntl[Controls::nbparts] : comm.size();
    }
    
    abcd::partitionMatrix();

    LINFO << "> Time to partition the matrix: "
          << MPI_Wtime() - t << "s.";
    
    t = MPI_Wtime();
    
    abcd::scaling();

    LINFO << "> Time to scale the matrix: "
          << MPI_Wtime() - t << "s.";

    abcd::analyseFrame();
    LINFO << "> Total time to preprocess: " << MPI_Wtime() - tot << "s.";
    LINFO << "*----------------------------------*";

    // everything is alright, tell the others that we're done here
    mpi::broadcast(comm, err, 0);

    return 0;
}

/// Creates augmented systems and factorizes them
int abcd::factorizeAugmentedSystems()
{
    // Create the group of CG instances
    abcd::createInterComm();

    double t = MPI_Wtime();
    
    if(instance_type == 0) {
        // some data on augmented system
        mpi::broadcast(inter_comm, m_o, 0);
        mpi::broadcast(inter_comm, n_o, 0);
        mpi::broadcast(inter_comm, nz_o, 0);

        mpi::broadcast(inter_comm, icntl, 0);
        mpi::broadcast(inter_comm, dcntl, 0);

        if( icntl[Controls::aug_type] != 0 )
            mpi::broadcast(inter_comm, size_c, 0);
    }

    if(instance_type == 0) {
        abcd::distributeData();
        abcd::createInterconnections();
    }
    
    abcd::initializeCimmino();
    
    if(inter_comm.rank() == 0 && instance_type == 0)
        LINFO << "Launching MUMPS analysis";
    
    abcd::analyseAugmentedSystems(mumps);
    
    if(IRANK == 0){
        LINFO << "Initialization time : " << MPI_Wtime() - t;
    }
    if(inter_comm.rank() == 0 && instance_type == 0){
        LINFO << "Launching MUMPS factorization";
    }

    t = MPI_Wtime();
    abcd::factorizeAugmentedSystems(mumps);
    
    if(IRANK == 0){
        LINFO << "Factorization time : " << MPI_Wtime() - t;
    }

    return 0;
}

/// Runs either BCG or ABCD solve depending on what we want
int abcd::solveSystem()
{
    if (nrhs < 1) {
        info[Controls::status] = -9;
        throw std::runtime_error("Number of right-hand sides should be at least one (1)");
    }
    if(icntl[Controls::block_size] < 1) {
        info[Controls::status] = -10;
        throw std::runtime_error("Block size should be at least one (1)");
    } 
    
    if(instance_type == 0) inter_comm.barrier();
    if(inter_comm.rank() == 0 && instance_type == 0){
        LINFO << "Launching Solve";
    }
    if(instance_type == 0) {

        mpi::broadcast(inter_comm, icntl, 0);
        mpi::broadcast(inter_comm, dcntl, 0);

        if(size_c == 0 && inter_comm.rank() == 0 && icntl[Controls::aug_type] != 0){
            LINFO << "Size of S is 0, therefore launching bcg";
            icntl[Controls::aug_type] = 0;
        }

        // If we are running ABCD, we have to ensure that we have no block-size > 1
        if(icntl[Controls::aug_type] == 0 ||
#ifdef WIP
           icntl[Controls::aug_project] != 0 ||
#endif //WIP
           runSolveS){
            abcd::distributeRhs();
            abcd::bcg(B);
        } else{
            int temp_bs = icntl[Controls::block_size];
            icntl[Controls::block_size] = 1;
            abcd::distributeRhs();
            abcd::solveABCD(B);
            icntl[Controls::block_size] = temp_bs;
        }

        int job = -1;
        mpi::broadcast(intra_comm, job, 0);

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
}
    
///  The gateway function that launches all other options
///
/// \param job The job id
/// \return Status code
/// @TODO Replace the number of jobs by identifiers easy to remember
int abcd::operator()(int job)
{
    info[Controls::status] = 0;

    LDEBUG3 << "MPI-Process " << comm.rank() << " called job = " << job;
    
    if ( job == 1 && last_called_job != -1 ) {
        info[Controls::status] = -2;
        throw std::runtime_error("Did you forget to call job = -1? ");
    }

    if ( job == 2 && job <= 3 && last_called_job != 1) {
        info[Controls::status] = -2;
        throw std::runtime_error("Did you forget to call job = 1? ");
    }

    if ( job == 3 && last_called_job != 2 && last_called_job != 3 ) {
        info[Controls::status] = -2;
        throw std::runtime_error("Did you forget to call job = 2? ");
    }

    if ( job < last_called_job && job != 3) {
        info[Controls::status] = -2;
        throw std::runtime_error("You cannot go back in time");
    }

    if ( job == last_called_job && job != 3) {
        info[Controls::status] = -2;
        throw std::runtime_error("You cannot re-run the same job unless you want to re-run job_id = 3 ");
    }
    
    switch(job) {

    case -1:
        // in case the user decided to change the log_file name or
        // decided to disable the logging
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

        initializeMatrix();

        // if everything went alright, remember the job
        last_called_job = job;

        break;

    case 1:
        preprocessMatrix();

        // if everything went alright, remember the job
        last_called_job = job;

        break;

    case 2:
        factorizeAugmentedSystems();

        // if everything went alright, remember the job
        last_called_job = job;

        break;

    case 3:
        solveSystem();

        // if everything went alright, remember the job
        last_called_job = job;

        break;

    case 4:
        // this ensures that we get consistent results
        operator()(1);
        operator()(2);

        // if everything went alright, remember the job
        last_called_job = 2;

        break;

    case 5:
        operator()(2);
        operator()(3);

        // if everything went alright, remember the job
        last_called_job = 3;

        break;

    case 6:
        operator()(1);
        operator()(2);
        operator()(3);

        // if everything went alright, remember the job
        last_called_job = 3;

        break;

    default:
        // Wrong job id
        info[Controls::status] = -1;
        throw std::runtime_error("Wrong job id.");
    }
    return 0;
}
