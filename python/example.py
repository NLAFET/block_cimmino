#!/usr/bin/env python2

import mpi4py.MPI as MPI
import numpy as np
from scipy.io import mmread
from scipy.linalg import solve
import abcdpy
from abcdpy import abcd
from time import time

solver = abcd()
A = mmread("/home/mzenadi/projects/abcd/example/e05r0500.mtx").tocoo()
b = mmread("/home/mzenadi/projects/abcd/example/e05r0500_rhs.mtx")

t = time()

solver.set_matrix(A.shape[0], A.shape[1], A.nnz, A.row, A.col, A.data)
solver.set_rhs(b, b.shape[1])

# augment the matrix using A_ij/-A_ji style
solver.icntl[abcdpy.icontrols.aug_type] = 2

#be a little bit verbose
solver.icntl[abcdpy.icontrols.verbose_level] = 1

solver.run(-1)
solver.run(6)

print "Backward error : {0:.2e}".format(solver.dinfo[abcdpy.dinfo.backward])
print "Took {0:.3f} seconds to solve".format(time() - t)

# get the solution vector
solution = solver.get_sol()

# plot S in the case of augmentation
if solver.icntl[abcdpy.icontrols.aug_type] != 0:
    sr = solver.get_s()
    import ipdb; ipdb.set_trace()
