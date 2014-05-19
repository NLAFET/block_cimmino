.. ABCD Solver documentation master file, created by
   sphinx-quickstart on Mon May 19 18:01:08 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**********************************************
The Augmented Block Cimmino Distributed Solver
**********************************************

.. toctree::
   :maxdepth: 2

=============
The Algorithm
=============

The package `ABCD Solver` is a distributed hybrid (iterative/direct)
solver for sparse linear systems :math:`Ax = b` where :math:`A` is a double
precision matrix with. `ABCD Solver` uses two methods to solve the linear system:

- *Regular Block Cimmino*: A block-projection technique that iterates
   to solve the linear system. During the iterations it solves a set
   of small problems (augmented systems built using the partitions of
   the original system).
- *Augmented Block Cimmino*: A pseudo-direct technique that augments
   the original system and through a succession of direct solves finds
   the solution.

The regular block Cimmino
-------------------------

The block Cimmino method is an iterative method that uses block-row
projections. To solve the :math:`Ax = b`, where :math:`A` is an
:math:`m\times n` sparse matrix, :math:`x` is an :math:`n`-vector and
:math:`b` is an :math:`m`-vector, we subdivide the system into strips of
rows as follows:

.. math::
   :nowrap:

   \begin{eqnarray*}
    \left(
      \begin{array}{c}
        A_1 \\ A_2 \\ \vdots \\ A_p
      \end{array}
    \right)
    ~x &=&
    \left(
      \begin{array}{c}
        b_1 \\ b_2 \\ \vdots \\ b_p
      \end{array}
    \right)
   \end{eqnarray*}


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

