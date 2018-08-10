// Copyright Institut National Polytechnique de Toulouse (2014)
// Contributor(s) :
// M. Zenadi <mzenadi@enseeiht.fr>
// D. Ruiz <ruiz@enseeiht.fr>
// R. Guivarch <guivarch@enseeiht.fr>

// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use,
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html"

// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability.

// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.

// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.

/*!
 * \file abcd_c.h
 * \brief Interface of the ABCD solver for the C programming language
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#ifndef _ABCD_C_H_
#define _ABCD_C_H_

enum controls {
   abcd_nbparts            ,
   abcd_part_type          ,
   abcd_part_guess         ,
   abcd_scaling            ,
   abcd_itmax              ,
   abcd_block_size         ,
   abcd_verbose_level      ,
   abcd_aug_type           ,
   abcd_aug_blocking       ,

   abcd_part_imbalance     ,
   abcd_threshold          ,

   abcd_status             ,
   abcd_nb_iter            ,
   abcd_residual           ,
   abcd_forward_error      ,
   abcd_backward           ,
   abcd_scaled_residual
};


struct abcd_solver
{
  int m; ///< row number
  int n; ///< column number
  int nz; ///< number of nnz in the lower-triangular part
  int *icntl;
  int *info;
  double *dcntl;
  double *dinfo;
  int sym;
  char *write_problem;

  int *irn;
  int *jcn;
  double *val;
  int start_index;

  int nrhs;
  double *rhs;
  double *sol;
};

typedef struct abcd_solver abcd_c;

#ifndef __cplusplus
struct abcd_solver* new_solver();
    void call_solver(struct abcd_solver *solver, int job_id);
#else
extern "C" {
    extern struct abcd_solver* new_solver();
    extern void call_solver(struct abcd_solver *solver, int job_id);
}
#endif

#endif // _ABCD_C_H_
