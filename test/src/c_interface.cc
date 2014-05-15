#include "abcd.h"
#include "abcd_c.h"
#include <string>

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
  //obj->write_problem = string(solver->write_problem);
  obj -> m = solver -> m;
  obj -> n = solver -> n;
  obj -> nz = solver -> nz;
  obj -> sym = (solver -> sym != 0);
  obj -> irn = solver -> irn;
  obj -> jcn = solver -> jcn;
  obj -> val = solver -> val;
  obj -> rhs = solver -> rhs;
  cout << obj-> sym << endl;
  cout << obj-> m << endl;
  cout << obj-> n << endl;
  cout << obj-> nz << endl;
  for (size_t i = 0; i < obj -> nz; i++) {
    cout << obj-> irn[i] << "  " << obj-> jcn[i] << "  " << obj -> val[i]<< endl;
  }  

  // TRAVAIL
  obj -> operator()(job_id);
  //cout << obj->write_problem << endl;
  //obj -> write_problem =  string("aaa");
  // CPP -> C
  //solver -> write_problem = (char *)(obj -> write_problem).c_str();
}

