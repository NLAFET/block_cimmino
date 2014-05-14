#include <abcd.h>
#include <algorithm>
#include <exception>

#include <ctime>
using namespace std;

_INITIALIZE_EASYLOGGINGPP

/// Set the defaults in the constructor
abcd::abcd()
{
    last_called_job = 0;
    
    nrhs = 1;
    start_index = 0;
    use_xk = false;
    use_xf = false;
    rhs = NULL;
    size_c = 0;
    verbose = false;
    runSolveS = false;

    irn = nullptr;
    jcn = nullptr;
    val = nullptr;

    icntl.assign(20, 0);
    dcntl.assign(20, 0);
    info.assign(2, 0);
    dinfo.assign(5, 0);

    icntl[Controls::aug_blocking] = 256;

    icntl[Controls::nbparts] = 4;
    icntl[Controls::part_type] = 2;
    dcntl[Controls::part_imbalance] = 0.5;
    icntl[Controls::part_guess] = 1;
    icntl[Controls::scaling] = 2;
    icntl[Controls::itmax] = 1000;
    icntl[Controls::block_size] = 1;
    dcntl[Controls::threshold] = 1e-12;

    icntl[Controls::verbose_level] = 0;
    info[Controls::status] = 0;

    mpi::communicator world;
    comm = world;


    // Prepare for Logging!
    if (comm.rank() == 0 && log_output == "") {
        time_t t = time(0);
        struct tm * now = localtime(&t);

        char timeString[15];  // space for "HH_MM_SS_dd_mm\0"

        strftime(timeString, sizeof(timeString), "%H_%M_%S_%d_%m", now);

        std::string tt(timeString);
        log_output = "abcd_" + tt + ".log";
    }
    configure_logger(log_output);
}

abcd::~abcd()
{
}

/// Creates the internal matrix from user's data
int abcd::initializeMatrix()
{
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
    LINFO << "M  = " << m;
    LINFO << "N  = " << n;
    LINFO << "NZ = " << nz;

    if(m <= 0 || n <= 0 || nz <= 0){
        info[Controls::status] = -2;
        LOG_IF(m <= 0, ERROR) << "m is negative or zero";
        LOG_IF(n <= 0, ERROR) << "n is negative or zero";
        LOG_IF(nz<= 0, ERROR) << "nz is negative or zero";
        throw std::range_error("Errornous information about the matrix");
    }

    double t = MPI_Wtime();
    if(sym) {
        //over estimate nz
        int     *t_irn = new int[2*nz];
        int     *t_jcn = new int[2*nz];
        double  *t_val = new double[2*nz];
        int     t_nz = 0;
        for(int k = 0; k < nz; ++k) {
            irn[k]--; jcn[k]--;

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

        LINFO << "Local matrix initialized in " << MPI_Wtime() - t << "s.";
        delete[] t_irn;
        delete[] t_jcn;
        delete[] t_val;
    } else {
        t = MPI_Wtime();
        for(int i=0; i<nz; i++){
            irn[i]--;
            jcn[i]--;
        }
        Coord_Mat_double t_A;
        t_A = Coord_Mat_double(m, n, nz, val, irn, jcn, MV_Matrix_::ref);
        A = CompRow_Mat_double(t_A);
        LINFO << "Local matrix initialized in " << MPI_Wtime() - t << "s.";
    }
    
    n_o = n;
    m_o = m;
    nz_o = nz;
        
    LINFO << "Matrix initialization done";
    return 0; 
}

/// Scales, partitions and analyses the structure of partitions
int abcd::preprocessMatrix()
{

    if(comm.rank() != 0) return 0;

    nbparts = icntl[Controls::nbparts];
    if (parallel_cg == 0) {
        parallel_cg = nbparts < comm.size() ? nbparts : comm.size();
    }
    
    abcd::partitionMatrix();
    
    double timeToPreprocess = MPI_Wtime();
    abcd::scaling();
    abcd::analyseFrame();
    LINFO << "Time for preprocess : " << MPI_Wtime() - timeToPreprocess;
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
        if(icntl[Controls::aug_type] == 0 || icntl[Controls::aug_project] != 0 || runSolveS){
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
int abcd::operator()(int job)
{
    LDEBUG3 << "MPI-Process " << comm.rank() << " called job = " << job;
    
    if ( (job == 1 || job == 5 || job == 6) && last_called_job != -1 )
        throw std::runtime_error("Did you forget to call job = -1? ");
    if ( job == 2 && last_called_job != 1)
        throw std::runtime_error("Did you forget to call job = 1? ");
    if ( job == 3 && last_called_job != 2)
        throw std::runtime_error("Did you forget to call job = 2? ");
    
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
        break;

    case 1:
        preprocessMatrix();
        break;

    case 2:
        factorizeAugmentedSystems();
        break;

    case 3:
        solveSystem();
        break;

    case 5:
        preprocessMatrix();
        factorizeAugmentedSystems();
        break;

    case 6:
        preprocessMatrix();
        factorizeAugmentedSystems();
        solveSystem();
        break;

    default:
        // Wrong job id
        throw std::runtime_error("Wrong job id.");
        return -1;
    }

    // if everything went alright, remember the job
    last_called_job = job;
    return 0;
}
