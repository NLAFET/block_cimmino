struct abcd_solver
{
  int *icntl;
  int *dcntl;
};

#ifndef __cplusplus
struct abcd_solver* new_solver();
#else
extern "C" {
    extern struct abcd_solver* new_solver();
}
#endif
