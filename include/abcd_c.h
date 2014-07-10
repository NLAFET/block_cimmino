struct abcd_solver
{
  int m; ///< row number
  int n; ///< column number
  int nz; ///< number of nnz in the lower-triangular part
  int *icntl;
  int *dcntl;
  int sym;
  
  char *write_problem;

  int *irn;
  int *jcn;
  double *val;

  double *rhs;
};

#ifndef __cplusplus
struct abcd_solver* new_solver();
    void call_solver(struct abcd_solver *solver, int job_id);
#else
extern "C" {
    extern struct abcd_solver* new_solver();
    extern void call_solver(struct abcd_solver *solver, int job_id);
}
#endif
