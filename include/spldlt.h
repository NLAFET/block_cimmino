/*!
 * \file spldlt.h
 * \brief Implementation of interface to SpLDLT solver
 * \author S. Cayrols
 * \version 1.0
 * \date Oct 18th, 2018
 */

#ifndef _SpLDLT_HXX_
#define _SpLDLT_HXX_

//TODO uncomment when done in SpLDLT
#include "spldlt_iface.h" 
#include "spral_matrix_util.h"

struct SPLDLT : spldlt_data_t {
//struct SPLDLT {
public:
    bool initialized;
    int job;
    int n;        //#columns of the augmented system
    int nz;       //#non-zeros of the augmented system
    //Not figured out storage type of the augmented system yet!
    int *row;     //rows of the augmented system
    long *ptr;     //Columns of the augmented system
    double *val;    //Values of the augmented system
    int nrhs;     //#rhs
    int lrhs;     //Guessing leading dimension
    double *rhs;  //rhs
    int ncpu;
    int ngpu;

    SPLDLT() : initialized(false) {}
  //void operator()(int job_id) {
  //    this->job = job_id;
  //  //dmumps_c(this); //Replace if allowed by SpLDLT, otherwise to remove
  //}

};

#endif //_SpLDLT_HXX_
