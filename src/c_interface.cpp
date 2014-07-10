#include "abcd.h"
#include "abcd_c.h"
#include <string>
#include <algorithm>

abcd *obj;

#include <iostream>
using namespace std;
/** description breve
 *
 * detailleeeee
 * @return a new solver
 */
struct abcd_solver* new_solver()
{
    struct abcd_solver *solver = (struct abcd_solver*) malloc(sizeof(*solver));

    obj = new abcd();
    //solver -> icntl = &obj -> icntl[0];
    return solver;
}
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
    obj->write_problem = string(solver->write_problem);

    for(size_t i = 0; i < obj -> icntl.size(); i++){
        obj -> icntl[i] = solver -> icntl[i];
    }

    // call the solver
    obj -> operator()(job_id);

    // CPP -> C
    //solver -> write_problem = (char *)(obj -> write_problem).c_str();
}


void free_solver(struct abcd_solver *solver)
{
    delete obj;

    free(solver);
}
