=======================
Using the solver
=======================

The solver is in the form of a class named ``abcd`` and the user has to
instantiate it on each MPI-process to be involved.
In the following, we refer to the members of the class by ``member``
and to the methods by ``method()``. Arrays will have ``[]`` appended
to them, if we specify a size then the array is pre-allocated at
construction, otherwise it is either allocated by the user (such as
the linear system entries) or by the solver once it's generated (such
as the solution vector). The user can refer to :ref:`section_controls`
for more details.

Instantiating and calling the solver
------------------------------------

To use the solver, the user has to instantiate the class `abcd`.
During the construction of the instance, the default parameters are
initialized. The user can then call the object as a function
(*functor*) with the job number as an argument.

.. code-block:: cpp

    abcd obj; // instantiating the class
    obj(job_id); // call the solver with a job identifier job_id

For ``C``, the solver is a structure called ``struct abcd``:

.. code-block:: c

    structure abcd_solver *obj = new_solver(); // create a new solver
    call_solver(obj, job_id);

``job_id`` defines which operation the solver has to run. It can have
values **-1**, and **1** through **6**. The order in
which these jobs have to be called is described in :ref:`Job
dependencies <job_flow>`.

* **-1**, initializes the internal matrix used by the solver. Prior to this call, the user must provide:

  - The information about the matrix ``m``, ``n``, ``nz``,
    ``sym``, ``irn[]``, ``jcn[]``, ``val[]`` have to be initialized
    before the call. See [Input matrix and right-hand side] for more detail.
  - After the call, the arrays ``irn[]``, ``jcn[]``, ``val[]`` are no longer used by the solver.
* **1**, performs the preprocessing. During this call, the solver
  scales the matrix, partition it, and if required by the user
  performs the augmentation of the matrix. Prior to this call, the
  user must provide:

  - The number of partitions to create (see ``icntl[nbparts]``) or ask the solver to guess the appropriate number of partitions (see ``icntl[part_guess]``)
  - The type of scaling to perform (see ``icntl[scaling]``)
  - The type of augmentation to perform (see ``icntl[aug_type]``)
* **2**, creates the augmented systems, analyses them, creates the mapping between the different mpi-processes and factorizes the augmented systems.
* **3**, performs the solution step, the right-hand sides and their number are required prior to this call.

  - The right-hand sides have to be given through the array ``rhs[]`` and their number in ``nrhs``.
  - The block-size to be used during the bloc-CG acceleration. Its value is used only during the regular block Cimmino solve, and by default its value is 1.
  - The solution is centralized (on the master) in the array ``sol[]``.
* **4**, regroups the call to the phases 1 and 2. 
* **5**, regroups the call to the phases 2 and 3. 
* **6**, regroups the call to the phases 1, 2 and 3.

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

Input matrix and right-hand side
--------------------------------

The current version of the ABCD Solver accepts only real, centralized, linear systems. The definition of the linear system uses 7 members:

- ``m`` (type: ``int``), the number of rows.
- ``n`` (type: ``int``), the number of columns. 
- ``nz`` (type: ``int``), the number of entries.
- ``sym`` (type: ``bool``), the symmetry of the matrix. If the matrix is symmetric, the matrix must be given in a lower-triangular form.
- ``irn`` (type: ``int *``), the row indices. 
- ``jcn`` (type: ``int *``), the column indices.
- ``val`` (type: ``double *``), the matrix entries.
- ``rhs`` (type: ``double *``), the right-hand sides.
- ``nrhs`` (type: ``int``), the number of right-hand sides (default value is 1)

If either of the row and column indices start with **0** the arrays
are supposed to be zero based (`C` arrays indexation), otherwise, if
they start with **1** the arrays are supposed to be one based (`Fortran`
arrays indexation). If however, none starts with **0** or **1** then there
is either an empty row or an empty column and the solver raises an exception.

In ``C++``:

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

In ``C``:

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
``defaults.h``. To access them, the user can use the namespace
``Controls``, eg. ``Controls::scaling`` has a value of ``5`` and is used
with ``icntl`` to handle the scaling of the linear system.

The integer control array
#########################

- `abcd::icntl[Controls::nbparts]` or `abcd::icntl[1]` defines the number of partitions in our linear system, can be from `1` to `m` (the number of rows in the matrix)

.. code-block:: cpp

    // we have 8 partitions
    obj.icntl[Controls::nbparts] = 8;

- `abcd::icntl[Controls::part_type]` or `abcd::icntl[2]` defines the partitioning type. It can have the values:

    + `1`, manual partitioning, the *nbparts* partitions can be provided into the STL vector `obj.nbrows[]`. Example:
~~~~~~~~~~~~~~~{.cpp}
    // use manual partitioning
    obj.icntl[Controls::part_type] = 1;
    // say that we want 20 rows per partition
    obj.nrows.assign(obj.icntl[Controls::nbparts], 20);

    // or 
    obj.nrows.resize(obj.icntl[Controls::nbparts]);
    obj.nrows[0] = 20;
    obj.nrows[1] = 20;
    //...
~~~~~~~~~~~~~~~
    + `2` (*default*), automatic uniform partitioning, creates *nbparts* partitions of similar size.
~~~~~~~~~~~~~~~{.cpp}
    // use patoh partitioning
    obj.icntl[Controls::part_type] = 2;
~~~~~~~~~~~~~~~
    + `3`, automatic hypergraph partitioning, creates *nbparts* partitions using the hypergraph partitioner `PaToH`. The imbalance between the partitions is handled using `obj.dcntl[Controls::part_imbalance]`. Example:
~~~~~~~~~~~~~~~{.cpp}
    // use patoh partitioning
    obj.icntl[Controls::part_type] = 3;
    // say that we want an imbalance of 0.3 between the partitions
    obj.dcntl[Controls::part_imbalance] = 0.3;
~~~~~~~~~~~~~~~

- `abcd::icntl[Controls::part_guess]` or `abcd::icntl[4]` asks the solver to guess the appropriate number of partitions and overrides the defined *nbparts*. 

    + `0` **default**, no guess
    + `1`, guess

- `abcd::icntl[Controls::scaling]` or `abcd::icntl[5]` defines the type of scaling to be used.

    + `0`, no scaling
    + `1`, infinity norm `MC77` based scaling
    + `2` **default**, combination of one norm and two norm `MC77` based scaling
- `abcd::icntl[Controls::itmax]` or `abcd::icntl[6]` defines the maximum number of iterations in block-CG acceleration, default is `1000`
- `abcd::icntl[Controls::block_size]` or `abcd::icntl[7]` defines the block-size to be used by the block-CG acceleration, default is `1` for classical CG acceleration
- `abcd::icntl[Controls::verbose_level]` or `abcd::icntl[8]` defines how verbose the solver has to be. 
- `abcd::icntl[Controls::aug_type]` or `abcd::icntl[10]` defines the augmentation type.

    + `0` **default**, no augmentation. This makes the solver run in
    **regular block Cimmino** mode.
    + `1`, makes the solver run in **Augmented Block Cimmino** mode
    with an augmentation of the matrix using the \f$C_{ij}/-I\f$
    technique. For numerical stability, this augmentation technique
    has to be used with a scaling.
    + `2`, makes the solver run in **Augmented Block Cimmino** mode
    with an augmentation of the matrix using the \f$A_{ij}/-A_{ji}\f$
    technique. This is the prefered augmentation technique.

- `abcd::icntl[Controls::aug_blocking]` or `abcd::icntl[11]` defines the blocking factor when building the auxiliary matrix $S$, default is `128`. 
- `abcd::icntl[Controls::aug_analysis]` or `abcd::icntl[12]`, when set to a value different than `0`, analyses the number of columns in the augmentation.
- `abcd::icntl[13]` to `abcd::icntl[16]` are for development and testing purposes only.

## The double precision control array ## {#subsection_double_controls}
- `abcd::dcntl[Controls::part_imbalance]` or `obj.dcntl[1]` defines the imbalance between the partitions when using `PaToH` (`abcd::icntl[Controls::part_imbalance] = 3`).
- `obj.dcntl[Controls::threshold]` or `abcd::dcntl[2]` defines the stopping threshold for the block-CG acceleration, default is `1e-12`.

# A usage example # {#section_introduction}

Combining the previous options, we expose a basic example that uses
the regular block Cimmino scheme, we comment the interesting parts and
explain how they fit together.  
Refer to [Calling a job], [The linear system], and [The Controls] for more details.

~~~~~~~~~~~~~~~{.cpp}
    #include "abcd.h"
    // use boost::mpi for simplicity, the user can use which ever he wants
    #include "mpi.h"
    #include <boost/mpi.hpp>

    int main(int argc, char* argv[]) 
    {
        mpi::environment env(argc, argv);
        // obtain the WORLD communicator, by default the solver uses it
        mpi::communicator world;

        // create one instance of the abcd solver per mpi-process
        abcd obj;

        if(world.rank() == 0) { // the master
            // we create a 5x5 matrix for a 1D mesh + three-point stencil 
            obj.sym = true; // the matrix is symmetric
            obj.m = 10; // number of rows
            obj.n = obj.m; // number of columns
            obj.nz = 2*obj.m - 1; // number of nnz in the lower-triangular part

            // allocate the arrays
            obj.irn = new int[obj.nz];
            obj.jcn = new int[obj.nz];
            obj.val = new double[obj.nz];

            // initialize the matrix
            // Notice that the matrix is stored in 1-based format
            size_t pos = 0;
            for (size_t i = 1; i < obj.m; i++) {
                // the diagonal
                obj.irn[pos] = i;
                obj.jcn[pos] = i;
                obj.val[pos] = 2.0;
                pos++;

                // the lower-triangular part
                obj.irn[pos] = i + 1;
                obj.jcn[pos] = i;
                obj.val[pos] = -1.0;
                pos++;
            }

            // the last diagonal element
            obj.irn[pos] = obj.m;
            obj.jcn[pos] = obj.m;
            obj.val[pos] = 2.0;

            pos++;

            // set the rhs
            obj.rhs = new double[obj.m];
            for (size_t i = 0; i < obj.m; i++) {
                obj.rhs[i] = ((double) i + 1)/obj.m;
            }

            // ask the solver to guess the number of partitions
            obj.icntl[Controls::part_guess] = 1;
        }

        try {
            // We call the solver directly using the object itself
            // (the abcd class is a functor)
            obj(-1); // initialize the object with defaults
            obj(6); // equivalent to running 1, 2 and 3 successively
            // the solution is stored in obj.sol
        } catch (runtime_error err) {
            // In case there is a critical error, we throw a runtime_error exception
            cout << "An error occured: " << err.what() << endl;
        }
        
        if(world.rank() == 0) { // the master
            delete[] obj.irn;
            delete[] obj.jcn;
            delete[] obj.val;
        }

      return 0;
    }
~~~~~~~~~~~~~~~
