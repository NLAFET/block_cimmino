/*!
 * \file mumps.h
 * \brief Implementation of interface to MUMPS solver
 * \author M. Zenadi
 * \version 1.0
 * \date Aug 15, 2012
 */

#ifndef _MUMPS_HXX_
#define _MUMPS_HXX_

#include "dmumps_c.h"

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
