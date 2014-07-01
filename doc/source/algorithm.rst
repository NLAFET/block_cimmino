=============
The Algorithm
=============

The package **Augmented Block Cimmino Distributed Solver** (`ABCD
Solver`) is a distributed hybrid (iterative/direct) solver for sparse
linear systems :math:`Ax = b` where :math:`A` is a double precision
matrix :math:`m \times n`, :math:`x` an :math:`n`-vector and :math:`b`
an :math:`m`-vector.  `ABCD Solver` uses two methods to solve the
linear system:

- *Regular Block Cimmino*: A block-projection technique that iterates
   to solve the linear system. During the iterations it solves a set
   of small problems (augmented systems built using the partitions of
   the original system).
- *Augmented Block Cimmino*: A pseudo-direct solver that augments the
   original system and constructs the solution directly out of
   independent solves using the augmented subsystems.

The regular block Cimmino
-------------------------

The block Cimmino method is an iterative method that uses block-row
projections. To solve :math:`Ax = b`, where :math:`A` is an
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

Let :math:`P_{\mathcal{R}(A_i^T)}` be the projector onto the range of
:math:`A_i^T` and :math:`{A_i}^+` be the Moore-Penrose pseudo-inverse of the
partition :math:`A_i`. The block Cimmino algorithm then computes a solution
iteratively from an initial estimate :math:`x^{(0)}`, viz:

.. math::
   :nowrap:

    \begin{eqnarray*}
        \begin{array}{ccl}
        u_{i}  & = & A_i^+ \left ( {b_i - A_i x^{(k)}} \right ) ~~~ i = 1, .... p \\
        x^{(k+1)}  & = & x^{(k)} + \omega \sum_{i=1}^p{u_i}
        \end{array}
    \end{eqnarray*}

where we can see the independence of the set of :math:`p` equations,
which is why the method is so attractive in a parallel environment.

With the above notations, the iteration equation is thus:

.. math::
   :nowrap:

    \begin{eqnarray*}
        x^{(k+1)} & = & x^{(k)} + \omega \sum_{i=1}^p{A_i^+ \left ( {b_i - A_i x^{(k)}} \right )} \\
          & = & \left( {I - \omega \sum_{i=1}^p{A_i^+ A_i}} \right) x^{(k)}
        + \omega \sum_{i=1}^p{A_i^+ b_i} \\
          & = & Q x^{(k)} + \omega \sum_{i=1}^p{A_i^+ b_i}.
          \label{something}
    \end{eqnarray*}

The iteration matrix for block Cimmino method is :math:`H = I - Q`,
which corresponds to a sum of projectors :math:`H = \omega
\sum_{i=1}^p{\mathcal{P}_{\mathcal{R}(A_i^T)}}`. It is thus symmetric
and positive definite and so we can solve

.. math::
    H x ~=~ \xi,

where :math:`\xi = \omega \sum_{i=1}^p{A_i^+ b_i}`
using conjugate gradient or block conjugate gradient methods.  As the
relaxation scalar :math:`\omega` appears on both sides of the
equation, we can set it to one.

At each step of the conjugate gradient algorithm we must solve for the
:math:`p` projections viz.

.. math::
   :nowrap:

    \begin{equation}
        A_i u_i ~=~ r_i, ~~~~  (r_i = {b_i - A_i x^{(k)}}),~~~ i = 1, .... p.
    \end{equation}

In our implementation we choose to solve these equations using the augmented system approach

.. math::
   :nowrap:

    \begin{eqnarray*}
        \left ( \begin{array}{cc} I & A_i^T \\ A_i & 0 \end{array} \right )
          \left ( \begin{array}{l} u_i \\ v_i \end{array} \right )
        &=&  \left ( \begin{array}{l} 0 \\ r_i \end{array} \right )
    \end{eqnarray*}

that we perform, at each iteration, using a direct method, giving
:math:`u_i = A_i^+ r_i` the projection needed for each partition
:math:`A_i`.  We use the multifrontal parallel solver :math:`MUMPS` to
do this.

Running our solver in the regular mode will go through the following steps:

- Partition the system into strips of rows (:math:`A_i` and :math:`b_i` for :math:`i = 1, \dots p`)
- Create the augmented systems
- Analyze and factorize the augmented systems using the direct solver :math:`MUMPS`
- Run a block conjugate gradient with an implicit iteration matrix
  :math:`H`, which resumes into independent augmented systems direct
  solves at each iteration.


The augmented block Cimmino
---------------------------

To understand the algorithm, suppose that we have a matrix :math:`A` with three partitions, described as follows:

.. math::
   :nowrap:
      
    \begin{equation}
        A =
        \left[
        \begin{array}{cccccc}
            A_{1,1} & A_{1,2} &&&&  A_{1,3}\\
            & A_{2,1} & A_{2,2} & A_{2,3} & \\
            &&& A_{3,2} & A_{3,3} &  A_{3,1}
        \end{array}
        \right].
    \end{equation}

Where :math:`A_{i,j}` is the sub-part of :math:`A_i`, the :math:`i`-th partition, that is interconnected algebraically to the partition :math:`A_j`, and vice versa.

The goal of the augmented block Cimmino algorithm is to make these
three partitions mutually orthogonal to each other, meaning that the
inner product of each couple of partitions is zero. We consider two
different ways to augment the matrix to obtain these zero matrix inner
products.

  * The first way to augment the matrix to make all the partitions mutually orthogonal to each other is obtained by putting the product :math:`C_{ij} = A_{ij}A_{ji}^T` on the right of the partition :math:`A_i` and adding :math:`-I` on the right of :math:`A_j` viz.

  .. math::
    \bar{A} =
    \left[
    \begin{array}{cccccc|ccc}
        A_{1,1} & A_{1,2} &         &          & A_{1,3} &         & C_{1,2}  & C_{1,3} &        \\
                & A_{2,1} & A_{2,2} & A_{2,3}  &         &         & -I       &         & C_{2,3}\\
                &         &         & A_{3,2}  & A_{3,3} & A_{3,1} &          & -I      & -I
    \end{array}\right].

    
  * The second way is to repeat the submatrices :math:`A_{ij}` and
    :math:`A_{ji}`, reversing the signs of one of them to obtain the
    augmented matrix :math:`\bar{A}` as in the following

  .. math::
    \bar{A} =
    \left[
    \begin{array}{cccccc|ccc}
        A_{1,1} & A_{1,2} &         &          & A_{1,3} &         & A_{1,2}  & A_{1,3} &        \\
                & A_{2,1} & A_{2,2} & A_{2,3}  &         &         & -A_{2,1} &         & A_{2,3}\\
                &         &         & A_{3,2}  & A_{3,3} & A_{3,1} &          & -A_{3,1}& -A_{3,2}
    \end{array}\right].

  This way :math:`\bar{A}_i\bar{A}_j^T` is zero for any pair :math:`i/j`, hence the new matrix has mutually orthogonal partitions.

Notice that we augment the matrix upper-down and shift the
augmentation at each step. This is done such that we do not create any new
interconnections between the new partitions. A simple check shows that
:math:`\bar{A}_i \bar{A}_j^T` is zero for any pair :math:`i/j`.
