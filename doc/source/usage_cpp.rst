=======================
Using the solver
=======================

Instantiating the solver
------------------------
To use the solver, the user has to instantiate the class ``abcd`` for
C++. In the case of C, the user creates a structure object using the
function ``new_solver()``.  During the construction of the instance,
the control parameters are initialized to their default value, see
:ref:`the controls description <section_controls>` for the list of the
control parameters and their default value.

For C++:

.. code-block:: cpp

    abcd obj; // instantiating the class

For C, the solver is a structure called ``struct abcd``:

.. code-block:: c

    structure abcd_solver *obj = new_solver(); // create a new solver

Input matrix and right-hand side
--------------------------------

The current version of the ABCD Solver accepts only real, centralized, linear systems. The definition of the linear system uses the following information:

.. doxygenclass:: abcd
    :project: abcd                  
    :members: m, n, nz, sym, irn, jcn, val, rhs, nrhs


If any of the row or column indices starts with **0** the arrays
are assumed to be zero based (`C` arrays indexation), otherwise, if
they start with **1** the arrays are supposed to be one based (`Fortran`
arrays indexation). If however, none starts with **0** or **1** then there
is either an empty row or an empty column and the solver raises an exception.

**Note** The solver does not check for empty rows or empty columns for the moment, and assumes that the given row indices are sorted.

To initialize the linear system in C++:

.. code-block:: cpp

    // Create an object for each mpi-process
    abcd obj;
    obj.n = 7;
    obj.m = 7;
    obj.nz = 15;
    obj.sym = false;
    if (world.rank() == 0) { // only the master is required to provide the matrix
        // allocate the arrays
        obj.irn = new int[obj.nz]
        // put the data in the arrays
        obj.irn[0] = 1;
        //..
    }

In C:

.. code-block:: cpp

    // Create an object for each mpi-process
    structure abcd_solver *obj = new_solver();

    obj->n = 7;
    obj->m = 7;
    obj->nz = 15;
    obj->sym = 0; 

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0) { // only the master is required to provide the matrix
        // allocate the arrays
        obj->irn = (int*) malloc(sizeof(int)*(obj->nz));
        // put the data in the arrays
        obj->irn[0] = 1;
        //..
    }

Calling the solver
------------------------------------
During the construction of the solver object, the default parameters are
initialized. The user can then call the object as a function
(*functor*) with the job number as an argument.

.. code-block:: cpp

    abcd obj; // instantiating the class
    obj(job_id); // call the solver with a job identifier job_id

For C, the solver is a structure called ``struct abcd``:

.. code-block:: c

    structure abcd_solver *obj = new_solver(); // create a new solver
    call_solver(obj, job_id);

``job_id`` defines which operation the solver has to run. It can have
values **-1**, and **1** through **6**. The order in
which these jobs have to be called is described in :ref:`Job
dependencies <job_flow>` figure.

.. doxygenclass:: abcd
    :project: abcd                  
    :members: operator()

.. _job_flow:

.. tikz:: Job dependencies.
   :libs: shapes, arrows, calc

    \tikzstyle{block} = [rectangle, draw, fill=blue!20, 
    text width=8.5em, text centered, minimum height=5em,
    node distance=12em]
    \tikzstyle{line} = [draw, -latex']
    \footnotesize

    \node [block] (init) {\texttt{Initialization (-1)}};
    \node [block, right of=init] (preprocess) {\texttt{Preprocessing (1)}};
    \node [block, right of=preprocess] (augsys) {\texttt{Analysis and Factorization (2)}};
    \node [block, right of=augsys] (solve) {\texttt{Solve (3)}};
    \node (retsol) at ($(augsys)!0.5!(solve)$) {};
    % Draw edges
    \path [line] (init) -- (preprocess);
    \path [line] (preprocess) -- (augsys);
    \path [line] (augsys) -- (solve);
    \path [line] (solve) -- +(0,-1.5) -| (retsol);

.. _section_controls:

The Controls
------------

Define the general behavior of the solver. They are split into two
arrays, `icntl` and `dcntl`. `icntl` is an *integer* array and defines
the options that control the specific parts of the solver, such as the
scaling, the type of algorithm to run and so on. `dcntl` is a *double
precision* array and defines some of the options required by the
algorithms we use such as the imbalance between the partition sizes
and the stopping criteria of the solver.

To access each of the control options we can either use the indices
``0, 1, ..`` or, preferably, use the *enums* defined in the header
``defaults.h`` (for ``C++``, in the case of ``C`` they are already
defined in ``abcd_c.h``). To access them, the user can use the
namespace ``Controls``, eg. ``Controls::scaling`` has a value of ``5``
and is used with ``icntl`` to handle the scaling of the linear system.
In the following we omit the ``Controls::`` part for simplicity.
Moreover, due to the similarities between the ``C++`` code and the
``C`` one, we provide only the ``C++`` snippets unless there is a
major difference.

The integer control array
#########################

.. doxygenenum:: icontrols
    :project: abcd

The double precision control array
##################################

.. doxygenenum:: dcontrols
    :project: abcd

* ``dcntl[part_imbalance]`` or ``obj.dcntl[1]`` defines the imbalance between the partitions when using ``PaToH`` (``icntl[part_imbalance] = 3``).
* ``obj.dcntl[threshold]`` or ``dcntl[2]`` defines the stopping threshold for the block-CG acceleration, default is ``1e-12``.

A usage example (C++)
---------------------

Combining the previous options, we expose a basic example that uses
the regular block Cimmino scheme, we comment the interesting parts and
explain how they fit together.  
Refer to [Calling a job], [The linear system], and [The Controls] for more details.

.. literalinclude:: ../../example/example.cpp
                    :language: cpp
