/*!
 * \file c_interface.cpp
 * \brief Implementation of the interface of ABCD for C
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include "abcd.h"
#include "abcd_c.h"
#include <string>
#include <algorithm>

abcd *obj;

#include <iostream>
using namespace std;

/*!
 *  \brief Interface to constructor of class ABCD
 *
 *  Interface to call the default constructor of the abcd class which sets
 *  default values for attributes.
 *
 */
struct abcd_solver* new_solver()
{
    struct abcd_solver *solver = (struct abcd_solver*) malloc(sizeof(*solver));

    obj = new abcd();

    solver->icntl = &obj->icntl[0];
    solver->dcntl = &obj->dcntl[0];
    solver->info  = &obj->info[0];
    solver->dinfo = &obj->dinfo[0];

    solver->irn         = obj->irn;
    solver->jcn         = obj->jcn;
    solver->val         = obj->val;
    solver->rhs         = obj->rhs;
    solver->nrhs        = obj->nrhs;
    solver->start_index = obj->start_index;
    solver->sym         = obj->sym;

    return solver;
}               /* -----  end of function new_solver  ----- */

/*!
 *  \brief Interface to the gateway function
 *
 *  Interface to the gateway function which launches all other options.
 *
 *  \param solver: structure of the ABCD solver
 *  \param job_id: the parameter choosing which phase of the solver to launch
 */
void call_solver(struct abcd_solver* solver, int job_id)
{
    // C -> CPP
    obj->m            = solver->m;
    obj->n            = solver->n;
    obj->nz           = solver->nz;
    obj->sym          = (solver->sym != 0);
    obj->irn          = solver->irn;
    obj->jcn          = solver->jcn;
    obj->val          = solver->val;
    obj->rhs          = solver->rhs;
    obj->nrhs         = solver->nrhs;
    obj->start_index  = solver->start_index;
    obj->parallel_cg  = solver->parallel_cg;

    // obj->write_problem = string(solver->write_problem);

    for(size_t i = 0; i < obj -> icntl.size(); i++){
        obj -> icntl[i] = solver -> icntl[i];
    }
    for(size_t i = 0; i < obj -> dcntl.size(); i++){
        obj -> dcntl[i] = solver -> dcntl[i];
    }

    // call the solver
    obj -> operator()(job_id);

    // CPP -> C
    //solver -> write_problem = (char *)(obj -> write_problem).c_str();
    solver->rhs = obj->rhs;
    solver->sol = obj->sol;
    solver->m   = obj->m;
    solver->n   = obj->n;
    solver->nz  = obj->nz;
    solver->sym = obj->sym ? 1 : 0;
}               /* -----  end of function call_solver  ----- */

/*!
 *  \brief Interface to the default Destructor of abcd class
 *
 *  Interface to the default destructor of the abcd class which deinitiliaze MUMPS.
 *
 */
void free_solver(struct abcd_solver *solver)
{
    delete obj;

    free(solver->icntl);
    free(solver->dcntl);
    free(solver->write_problem);
    free(solver);
}               /* -----  end of function free_solver  ----- */

/*!
 *  \brief Interface to the routine that parse the config_file
 *
 *  Interface to the routine that parse the config_file
 *
 */
int parse_configFile( struct abcd_solver  *solver,
                      char                *config_file,
                      char                **matrix_file,
                      char                **rhs_file){

  int err = 0;
  string matrix_file_str, rhs_file_str;

  err = obj -> parse_configFile(string(config_file),
      matrix_file_str, rhs_file_str);
  if(err) return err;

  /* Interface to save the problem */
  solver -> write_problem = (char*)malloc(obj -> write_problem.size());
  strncpy(solver -> write_problem,
          obj -> write_problem.c_str(),
          obj -> write_problem.size());
/*solver -> log_output = malloc(obj -> log_output.size());
  strncpy(solver -> log_output,
          obj -> log_output.c_str(),
          obj -> log_output.size());*/
  /* Interface matrix */
  *matrix_file = (char*)malloc(matrix_file_str.size());
  strncpy(*matrix_file,
          matrix_file_str.c_str(),
          matrix_file_str.size());
  /* Interface RHS */
  if(rhs_file_str.size()){
    *rhs_file = (char*)malloc(rhs_file_str.size());
    strncpy(*rhs_file,
            rhs_file_str.c_str(),
            rhs_file_str.size());
  }
/*solver -> start_file = malloc(obj -> start_file.size());
  strncpy(solver -> start_file,
          obj -> start_file.c_str(),
          obj -> start_file.size());*/

  for(size_t i = 0; i < obj -> icntl.size(); i++){
    solver -> icntl[i] = obj -> icntl[i];
  }
  for(size_t i = 0; i < obj -> dcntl.size(); i++){
    solver -> dcntl[i] = obj -> dcntl[i];
  }

  solver -> parallel_cg = obj -> parallel_cg;
  return 0;
}               /* -----  end of function parse_configFile  ----- */

/*!
 *  \brief Interface to the routine that loads the matrix
 *
 *  Interface to the routine that loads the matrix
 *
 */
int load_MM(struct abcd_solver *solver, char *matrix_file){
  int err = 0;

  err = obj -> load_MM(string(matrix_file));

  solver -> m           = obj -> m;
  solver -> n           = obj -> n;
  solver -> nz          = obj -> nz;
  solver -> start_index = obj -> start_index;
  solver -> irn         = obj -> irn;
  solver -> jcn         = obj -> jcn;
  solver -> val         = obj -> val;

  return err;
}

/*!
 *  \brief Interface to the routine that loads the RHS
 *
 *  Interface to the routine that loads the RHS
 *
 */
int load_RHS(struct abcd_solver *solver, char *rhs_file){
  int err = 0;

  err = obj -> load_RHS(string(rhs_file));
  
  solver -> nrhs  = obj -> nrhs;
  solver -> rhs   = obj -> rhs;

  return err; 
}
