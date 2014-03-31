#include "abcd.h"
#include "abcd_c.h"

abcd *obj;

#include <iostream>
using namespace std;
struct abcd_solver* new_solver()
{
  struct abcd_solver *solver = (struct abcd_solver *)malloc(sizeof(*solver));
  obj = new abcd();
  cout << obj->icntl[Controls::nbparts] << '\t' << Controls::nbparts << endl;
  solver->icntl = &obj->icntl[0];
  return solver;
}
