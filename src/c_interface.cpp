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
    solver->info = &obj->info[0];
    solver->dinfo = &obj->dinfo[0];

    solver->irn = obj->irn;
    solver->jcn = obj->jcn;
    solver->val = obj->val;
    solver->rhs = obj->rhs;
    solver->nrhs = obj->nrhs;
    solver->start_index = obj->start_index;

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
    obj->m   = solver->m;
    obj->n   = solver->n;
    obj->nz  = solver->nz;
    obj->sym = (solver->sym != 0);
    obj->irn = solver->irn;
    obj->jcn = solver->jcn;
    obj->val = solver->val;
    obj->rhs = solver->rhs;
    obj->nrhs = solver->nrhs;
    obj->start_index = solver->start_index;

    // obj->write_problem = string(solver->write_problem);

    for(size_t i = 0; i < obj -> icntl.size(); i++){
        obj -> icntl[i] = solver -> icntl[i];
    }

    // call the solver
    obj -> operator()(job_id);

    // CPP -> C
    //solver -> write_problem = (char *)(obj -> write_problem).c_str();
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
