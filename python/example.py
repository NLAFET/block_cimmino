#!/usr/bin/env python2

import mpi4py.MPI as MPI
import numpy as np
from scipy.io import mmread
from scipy.linalg import solve
from abcdpy import abcd

solver.= abcd()
A = mmread("/hsolver.e/mzenadi/projects/abcd/example/e05r0500.mtx").tocoo()
b = mmread("/hsolver.e/mzenadi/projects/abcd/example/e05r0500_rhs.mtx")

solver.set_matrix(A.shape[0], A.shape[1], A.nnz, A.row, A.col, A.data)
solver.set_rhs(b, b.shape[1])

solver.run(-1)
solver.run(6)

solution = o.get_sol()
