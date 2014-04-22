The Augmented Block Cimmino Distributed Solver {#mainpage}
=============

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

# Introduction to the methods #

## The regular block Cimmino ##

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

Let $P_{\mathcal{R}(A_i^T)}$ be the projector onto the range of
$A_i^T$ and ${A_i}^+$ be the Moore-Penrose pseudo-inverse of the
partition $A_i$. The block Cimmino algorithm then computes a solution
iteratively from an initial estimate $x^{(0)}$ according to:
\f{eqnarray*}
\begin{array}{ccl}
u_{i}  & = & A_i^+ \left ( {b_i - A_i x^{(k)}} \right ) ~~~ i = 1, .... p \\
x^{(k+1)}  & = & x^{(k)} + \omega \sum_{i=1}^p{u_i}
\end{array}
\f}
where we see the independence of the set of $p$ equations, which is why
the method is so attractive in a parallel environment.

With the above notations, the iteration equations are thus:
$$
    \begin{array}{ccl}
        x^{(k+1)} & = & x^{(k)} + \omega \sum_{i=1}^p{A_i^+ \left ( {b_i - A_i x^{(k)}} \right )} \\
          & = & \left( {I - \omega \sum_{i=1}^p{A_i^+ A_i}} \right) x^{(k)}
        + \omega \sum_{i=1}^p{A_i^+ b_i} \\
          & = & Q x^{(k)} + \omega \sum_{i=1}^p{A_i^+ b_i}.
          \label{something}
    \end{array}
$$

The iteration matrix for block Cimmino, $H = I - Q$, is then 
a sum of projectors $H = \omega \sum_{i=1}^p{\mathcal{P}_{\mathcal{R}(A_i^T)}}$.  
It is thus symmetric and positive definite
and so we can solve 
$$
    H x ~=~ \xi, 
$$
where $\xi = \omega \sum_{i=1}^p{A_i^+ b_i}$
using conjugate gradient or block conjugate gradient methods.  As $\omega$ appears on both sides of the equation, we can set it to one.

At each step of the conjugate gradient algorithm we must solve for the 
$p$ projections viz.
$$
    A_i u_i ~=~ r_i, ~~~~  (r_i = {b_i - A_i x^{(k)}}),~~~ i = 1, .... p.
$$

In our approach we choose to solve these equations using the augmented system
\f{eqnarray*}
    \left ( \begin{array}{cc} I & A_i^T \\ A_i & 0 \end{array} \right )
      \left ( \begin{array}{l} u_i \\ v_i \end{array} \right )
    &=&  \left ( \begin{array}{l} 0 \\ r_i \end{array} \right )
\f}
that we will solve, at each iteration, using a direct method and gives $u_i = A_i^+ r_i$ the projection we need for the partition $A_i$.

 We use the multifrontal parallel solver (\texttt{MUMPS})
\cite{adek:01} to do this. The main other techniques for solving
equation (\ref{eqn:projections}) are using normal equations or a QR
factorization. The former has numerical and storage issues while the
latter lacks a good distributed solver. We avoid both problems with
our approach.

Running our solver in the regular mode will go through the following steps:

- Partition the system into strips of rows ($A_i$ and $b_i$ for $i = 1, \dots p$)
- Create the augmented systems
- Analyze and factorize the augmented systems using the direct solver `MUMPS`
- Run a block conjugate gradient with an implicit matrix $H$, and at each iteration compute the matrix-vector product as a sum of projections. These projects being a set of solves using the direct solver.

## The augmented block Cimmino ##

To understand the algorith, suppose that we have a matrix $A$ with three partitions, described as follow:
\f[
    A = 
    \begin{bmatrix}
        A_{1,1} & A_{1,2} &&&&  A_{1,3}\\
        & A_{2,1} & A_{2,2} & A_{2,3} && \\
        &&& A_{3,2} & A_{3,3} &  A_{3,1}\\
    \end{bmatrix}
\f]

Where $A_{i,j}$ is the sub-part of $A_i$, the $i$-th partition, that is interconnected algebraically to the partition $A_j$, and vice versa.

The goal of the augmented block Cimmino algorithm is to make these
three partitions mutually orthogonal to each other, meaning that the
product of each couple of partitions is zero. We consider two
different ways to augment the matrix to obtain these zero matrix products. 

- One can repeat the submatrices $A_{ij}$ and $A_{ji}$, reversing the signs of one of them to obtain the augmented matrix $\bar{A}$ as in the following
$$
    \bar{A} = 
    \begin{bmatrix}
        A_{1,1} & A_{1,2} &&&&  A_{1,3}  &\vline& A_{1,2}  & A_{1,3}\\
        & A_{2,1} & A_{2,2} & A_{2,3}  &&&\vline& -A_{2,1} && A_{2,3}\\
        &&& A_{3,2} & A_{3,3} &  A_{3,1} &\vline&& -A_{3,1} & -A_{3,2}\\
    \end{bmatrix}.
$$
Notice that we augment the matrix upper-down and shift the
augmentation at each step. This way, we do not create any new
interconnections between the new partitions. A simple check shows that
$\bar{A}_i \bar{A}_j^T$ is zero for any couple $i/j$.


