/*!
 * \file mumps.h
 * \brief Implementation of interface to MUMPS solver
 * \author M. Zenadi
 * \version 1.0
 * \date Aug 15, 2012
 */

#ifndef _MUMPS_HXX_
#define _MUMPS_HXX_

#ifdef USE_MUMPS
  #include "dmumps_c.h"
#else
  #include <stdio.h>
  #include <stdlib.h>
  typedef struct {
    int sym;
    int par;
    int job;
    int comm_fortran;
    int n;
    int nz;
    int *irn;
    int *jcn;
    double *a;
    int nz_loc;
    int *irn_loc;
    int *jcn_loc;
    double *a_loc;
    int *perm_in;
    int *sym_perm;
    double *rhs;
    double *rhs_sparse;
    double *sol_loc;
    int *irhs_sparse;
    int *irhs_ptr;
    int *isol_loc;
    int nrhs;
    int lrhs;
    int nz_rhs;
    int lsol_loc;
    int icntl[40];
    int cntl[15];
    int info[40];
    int infog[40];
    int rinfo[40];
    int rinfog[40];
    char write_problem[256];
  } DMUMPS_STRUC_C;

  void dmumps_c(DMUMPS_STRUC_C * mu){
    printf("************\nMUMPS is not linked : "
        "please define the flag USE_MUMPS and recompile.\n");
    abort();
  }
#endif

struct MUMPS : DMUMPS_STRUC_C {
public:
    bool initialized;
    MUMPS() : initialized(false) {}
    void operator()(int job_id) {
        this->job = job_id;
        dmumps_c(this);
    }

    inline void setIcntl(int i, int v) { this->icntl[ i - 1 ] = v ; }
    inline void setCntl(int i, double v) { this->cntl[ i - 1 ] = v ; }

    inline const int getIcntl(int i) const { return this->icntl[ i - 1 ]; }
    inline const double getCntl(int i) const { return this->cntl[ i - 1 ]; }

    inline const int getInfo(int i) const { return this->info[ i - 1 ]; }
    inline const int getInfoG(int i) const { return this->infog[ i - 1 ]; }
    inline const double getRinfo(int i) const { return this->rinfo[ i - 1 ]; }
    inline const double getRinfoG(int i) const { return this->rinfog[ i - 1 ]; }
};

#endif //_MUMPS_HXX_
