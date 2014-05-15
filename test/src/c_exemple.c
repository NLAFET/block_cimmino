#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "abcd_c.h"

typedef struct abcd_solver abcd; // pkoi je peux pas le déclarer dans abcd.h ?

void init_2d_lap(abcd *o, int mesh_size);
void init_2d_lap2(int m, int n, int nz, int *irn, int *jcn, double *val, int mesh_size);

int main(int argc, char* argv[]) {

  int myrank;
  abcd *obj;

  int i;

  /* Initialize MPI */

  MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  obj = new_solver();

  if(myrank == 0) { // the master
    init_2d_lap(obj, 10);

    // set the rhs
    obj -> rhs = (double *) malloc(sizeof(double)*(obj -> m));
    for (i = 0; i < (obj -> m); i++) {
      obj -> rhs[i] = ((double) i + 1)/obj -> m;
    }
  }
  call_solver(obj, -1);
  call_solver(obj, 5);
 
  MPI_Finalize();

  return 0;
}

void init_2d_lap(abcd *obj, int mesh_size) {

    obj -> m = mesh_size*mesh_size; // number of rows
    obj -> n = obj -> m; // number of columns
    obj -> nz = 3*obj -> m - 2*mesh_size; // number of nnz in the lower-triangular part
    obj -> sym = 1;

    // allocate the arrays
    obj -> irn = (int*) malloc(sizeof(int)*(obj -> nz));
    obj -> jcn = (int*) malloc(sizeof(int)*(obj -> nz));
    obj -> val = (double*) malloc(sizeof(double)*(obj -> nz));

    init_2d_lap2(obj -> m, obj -> n, obj -> nz,
        obj -> irn, obj -> jcn, obj -> val, mesh_size);
    
    return;

}

void init_2d_lap2(int m, int n, int nz, 
    int irn[], int jcn[], double val[], int mesh_size) {

  int pos;
  int i;

  pos = 0;
  for (i = 1; i <= m; i++) {

    // the diagonal
    irn[pos] = i+1;
    jcn[pos] = i+1;
    val[pos] = 4.0;

    pos++;

    if (i == m) continue;
    // the lower-triangular part
    irn[pos] = i + 2;
    jcn[pos] = i + 1;
    val[pos] = -1.0;
    pos++;

    if (i > m - 2*mesh_size + 1) continue;
    irn[pos] = i + mesh_size + 1;
    jcn[pos] = i + 1;
    val[pos] = -1.0;
    pos++;
  }

  return;

}
