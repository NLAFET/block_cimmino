The Augmented Block Cimmino Distributed Solver {#mainpage}
==============

The package `ABCD Solver` is a distributed hybrid (iterative/direct)
solver for sparse linear systems \f$Ax = b\f$ where \f$A\f$ is a double
precision matrix with. `ABCD
Solver` uses two methods to solve the linear system:

- *Regular Block Cimmino*: A block-projection technique that iterates
   to solve the linear system. During the iterations it solves a set
   of small problems (augmented systems built using the partitions of
   the original system).
- *Augmented Block Cimmino*: A pseudo-direct technique that augments
   the original system and through a succession of direct solves finds
   the solution.

# The regular block Cimmino # {#section_cimmino}

The block Cimmino method is an iterative method that uses block-row
projections. To solve the \f$Ax = b\f$, where \f$A\f$ is an
\f$m\times n\f$ sparse matrix, \f$x\f$ is an \f$n\f$-vector and \f$b\f$ is an
\f$m\f$-vector, we subdivide the system into strips of rows as follows:

\f{eqnarray*}
  \left(
    \begin{array}{c}
      A_1 \\ A_2 \\ \vdots \\ A_p
    \end{array}
  \right)
  x &=&
  \left(
    \begin{array}{c}
      b_1 \\ b_2 \\ \vdots \\ b_p
    \end{array}
  \right).
\f}

Let \f$P_{\mathcal{R}(A_i^T)}\f$ be the projector onto the range of
\f$A_i^T\f$ and \f${A_i}^+\f$ be the Moore-Penrose pseudo-inverse of the
partition \f$A_i\f$. The block Cimmino algorithm then computes a solution
iteratively from an initial estimate \f$x^{(0)}\f$ according to:
\f{eqnarray*}
    \begin{array}{ccl}
    u_{i}  & = & A_i^+ \left ( {b_i - A_i x^{(k)}} \right ) ~~~ i = 1, .... p \\
    x^{(k+1)}  & = & x^{(k)} + \omega \sum_{i=1}^p{u_i}
    \end{array}
\f}
where we see the independence of the set of \f$p\f$ equations, which is why
the method is so attractive in a parallel environment.

With the above notations, the iteration equations are thus:
\f[
    \begin{array}{ccl}
        x^{(k+1)} & = & x^{(k)} + \omega \sum_{i=1}^p{A_i^+ \left ( {b_i - A_i x^{(k)}} \right )} \\
          & = & \left( {I - \omega \sum_{i=1}^p{A_i^+ A_i}} \right) x^{(k)}
        + \omega \sum_{i=1}^p{A_i^+ b_i} \\
          & = & Q x^{(k)} + \omega \sum_{i=1}^p{A_i^+ b_i}.
          \label{something}
    \end{array}
\f]

The iteration matrix for block Cimmino, \f$H = I - Q\f$, is then a sum
of projectors \f$H = \omega
\sum_{i=1}^p{\mathcal{P}_{\mathcal{R}(A_i^T)}}\f$. It is thus
symmetric and positive definite and so we can solve
\f[
    H x ~=~ \xi,
\f]
where \f$\xi = \omega \sum_{i=1}^p{A_i^+ b_i}\f$
using conjugate gradient or block conjugate gradient methods.  As \f$\omega\f$ appears on both sides of the equation, we can set it to one.

At each step of the conjugate gradient algorithm we must solve for the
\f$p\f$ projections viz.
\f{equation}
    A_i u_i ~=~ r_i, ~~~~  (r_i = {b_i - A_i x^{(k)}}),~~~ i = 1, .... p.
\f}

In our approach we choose to solve these equations using the augmented system
\f{eqnarray*}
    \left ( \begin{array}{cc} I & A_i^T \\ A_i & 0 \end{array} \right )
      \left ( \begin{array}{l} u_i \\ v_i \end{array} \right )
    &=&  \left ( \begin{array}{l} 0 \\ r_i \end{array} \right )
\f}
that we will solve, at each iteration, using a direct method and gives \f$u_i = A_i^+ r_i\f$ the projection we need for the partition \f$A_i\f$.
We use the multifrontal parallel solver `MUMPS` to do this.

Running our solver in the regular mode will go through the following steps:

- Partition the system into strips of rows (\f$A_i\f$ and \f$b_i\f$ for \f$i = 1, \dots p\f$)
- Create the augmented systems
- Analyze and factorize the augmented systems using the direct solver `MUMPS`
- Run a block conjugate gradient with an implicit matrix \f$H\f$, and at each iteration compute the matrix-vector product as a sum of projections. These projects being a set of solves using the direct solver.

# The augmented block Cimmino # {#section_abcd}

To understand the algorith, suppose that we have a matrix \f$A\f$ with three partitions, described as follow:

\f{equation}
    A =
    \left[
    \begin{array}{cccccc}
        A_{1,1} & A_{1,2} &&&&  A_{1,3}\\
        & A_{2,1} & A_{2,2} & A_{2,3} & \\
        &&& A_{3,2} & A_{3,3} &  A_{3,1}
    \end{array}
    \right].
\f}

Where \f$A_{i,j}\f$ is the sub-part of \f$A_i\f$, the \f$i\f$-th partition, that is interconnected algebraically to the partition \f$A_j\f$, and vice versa.

The goal of the augmented block Cimmino algorithm is to make these
three partitions mutually orthogonal to each other, meaning that the
product of each couple of partitions is zero. We consider two
different ways to augment the matrix to obtain these zero matrix products.

- The first way to augment the matrix to make all the partitions mutually orthogonal to each other is by putting the product \f$C_{ij} = A_{ij}A_{ji}^T\f$ on the right of the partition \f$A_i\f$ and put \f$-I\f$ on the right of \f$A_j\f$ viz.
\f{equation}
    \bar{A} =
    \left[
    \begin{array}{cccccc|ccc}
        A_{1,1} & A_{1,2} &         &          & A_{1,3} &         & C_{1,2}  & C_{1,3} &        \\
                & A_{2,1} & A_{2,2} & A_{2,3}  &         &         & -I       &         & C_{2,3}\\
                &         &         & A_{3,2}  & A_{3,3} & A_{3,1} &          & -I      & -I
    \end{array}\right].
\f}
This way \f$\bar{A}_i\bar{A}_j^T\f$ is zero for any pair \f$i/j\f$, hence the new matrix has mutually orthogonal partitions.

- We can repeat the submatrices \f$A_{ij}\f$ and \f$A_{ji}\f$, reversing the signs of one of them to obtain the augmented matrix \f$\bar{A}\f$ as in the following
\f[
    \bar{A} =
    \left[
    \begin{array}{cccccc|ccc}
        A_{1,1} & A_{1,2} &         &          & A_{1,3} &         & A_{1,2}  & A_{1,3} &        \\
                & A_{2,1} & A_{2,2} & A_{2,3}  &         &         & -A_{2,1} &         & A_{2,3}\\
                &         &         & A_{3,2}  & A_{3,3} & A_{3,1} &          & -A_{3,1}& -A_{3,2}
    \end{array}\right].
\f]
Notice that we augment the matrix upper-down and shift the
augmentation at each step. This way, we do not create any new
interconnections between the new partitions. A simple check shows that
\f$\bar{A}_i \bar{A}_j^T\f$ is zero for any pair \f$i/j\f$.


We shall call the first augmentation model the \f$C_{ij}/-I\f$
process, the second way as \f$A_{ij}/-A_{ji}\f$, and the augmentation
part as \f$C\f$. Implementation wise, they are enabled by setting
`abcd::icntl[Controls::aug_type]` to 1 for the former and 2 for the
latter.


