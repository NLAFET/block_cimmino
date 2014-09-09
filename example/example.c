#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "abcd_c.h"


void init_2d_lap(abcd_c *o, int mesh_size);
void init_2d_lap2(int m, int n, int nz, int *irn, int *jcn, double *val, int mesh_size);

int main(int argc, char* argv[]) {

  int myrank;
  abcd_c *obj;

  int i, j;

  /* Initialize MPI */

  MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  obj = new_solver();

  if(myrank == 0) { // the master
    init_2d_lap(obj, 200);

    obj->icntl[abcd_verbose_level] = 2;
    obj->icntl[abcd_scaling] = 2;
    obj->icntl[abcd_part_guess] = 1;

    // set the rhs
    printf("%d\n", obj->nrhs);
    
    obj->rhs = (double *) malloc(sizeof(double)*(obj->m * obj->nrhs));
    for (j = 0; j < obj->nrhs; ++j) 
        for (i = 0; i < obj->m; ++i) 
            obj->rhs[i + j * obj->m] = ((double) i * (j + 1))/obj->m;
  }
  call_solver(obj, -1);
  call_solver(obj, 6);
 
  MPI_Finalize();

  return 0;
}

void init_2d_lap(abcd_c *obj, int mesh_size) {

    obj->m = mesh_size*mesh_size; // number of rows
    obj->n = obj->m; // number of columns
    obj->nz = 3*obj->m - 2*mesh_size; // number of nnz in the lower-triangular part
    obj->sym = 1;
    obj->start_index = 1;

    // allocate the arrays
    obj->irn = (int*) malloc(sizeof(int)*(obj->nz));
    obj->jcn = (int*) malloc(sizeof(int)*(obj->nz));
    obj->val = (double*) malloc(sizeof(double)*(obj->nz));

    init_2d_lap2(obj->m, obj->n, obj->nz,
        obj->irn, obj->jcn, obj->val, mesh_size);
}

void init_2d_lap2(int m, int n, int nz, 
    int irn[], int jcn[], double val[], int mesh_size) {

  int pos;
  int i;

  pos = 0;
  for (i = 1; i <= m; i++) {

    // the diagonal
    irn[pos] = i;
    jcn[pos] = i;
    val[pos] = -4.0;

    pos++;

    if (i % mesh_size != 0) {
      // the lower-triangular part
      irn[pos] = i + 1;
      jcn[pos] = i ;
      val[pos] = 1.0;
      pos++;
    } 

    
    if (i + mesh_size > m) continue;
    irn[pos] = i + mesh_size ;
    jcn[pos] = i ;
    val[pos] = 1.0;
    pos++;
  }
}
