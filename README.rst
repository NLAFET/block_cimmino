==============================================
The Augmented Block Cimmino Distributed Solver
==============================================

**Note:** Check http://abcd.enseeiht.fr for more details.

Tested plateforms
-----------------

Working
=======

* Linux x86_64 with GNU 4.7 and 4.8  compilers,``MKL``, ``ACML``, ``OpenBLA`` reference blas and lapack.

Not Working
===========

* Fujitsu FX with Fujitsu compilers:

  - ``PaToH`` is not compatible (users have to request a compatible version from the authors)
  - Our Logging library is not compatible with Fujitsu compilers, should work with GNU compilers.

* Microsoft Windows:

  - ``MUMPS`` does not support Windows (there is an unofficial guide to compile it under Windows, but we 
do not provide any pre-compiled library for it)
  - ``PaToH`` is not compatible (users have to request a compatible version from the authors)

You can disable ``PaToH`` by running cmake with the option ``-DPATOH=OFF``. 

Not Tested
==========
* Mac OSX was not tested but should be fully compatible.    

Obtaining the source code
-------------------------

The ABCD Solver depends on a few external libraries: ``MUMPS``, ``Sparselib++ (custom)``, ``PaToH``, 
``lapack``, ``Boost::MPI`` and ``Boost::Serialization`` version 1.50 or higher.

* ``Sparselib++ (custom)``: a modified version of ``SparseLib++`` to
  suits our needs, is also distributed with our solver in the
  ``lib/sparselib`` directory. The library is already compiled,
  but you still can recompile it by running ``make all`` in
  ``lib/sparselib`` directory.
* ``MUMPS`` is mandatory. It is recommended to use the latest version 
  (currently 5.1.2) but any version ulterior to 5.1.0 is okay. MUMPS is available
  freely on demand on the MUMPS consortium website "mumps-solver.org".
  For versions older than 5.1.0, add the compilation flag "-DOLD_MUMPS" to 
  CMAKE_C_FLAGS/CMAKE_CXX_FLAGS in the CMakeLists.txt file in order to avoid issues 
  in the scaling part.
* ``PaToH``: Can be downloaded from the webpage of `Ümit V. Çatalyürek
  <http://bmi.osu.edu/~umit/software.html>`_ (URL available in the
  following script). The path ``libpatoh.a`` and the header `patho.h` must be respectively
  provided in LIB_DIRS and INC_DIRS of the abcdCmake configuration file.
* ``BLAS`` and ``LAPACK`` are both mandatory. We provide
  configurations to build the solver using ``ACML``, ``MKL`` and ``OpenBLAS``.
* ``BLACS`` and ``ScaLAPACK`` are required by ``MUMPS``, therefore
  they are needed when you link your software with the solver. We
  explicitly require them so that we can build the examples.
* ``Boost::MPI`` and ``Boost::Serialization`` are provided with the solver
  in the ``lib`` folder and are compiled at the same time as the abcd library.
  You can still provide an external Boost library by specifying ``BOOST_ROOT`` in
  the configuration file abcdCmake.in. Be careful that we do not guaranty compatibility
  with every version of Boost.
* ``Boost::MPI`` requires ``MPI`` and so does ``MUMPS``. You can
  install it either from source or through your distribution
  repositories. The solver was tested with versions 1.47, 1.49 and
  1.54. However, we recommend to use versions higher than 1.50.
	
The installation can be done by typing the following commands in your terminal

.. code-block:: bash

    # download the latest stable version
    # it will create a directory named abcd

    git clone https://bitbucket.org/apo_irit/abcd.git

    # download the appropriate version of patoh from
    # http://bmi.osu.edu/~umit/software.html

Now that everything is ready, we can compile the solver. To do so, we
need a configuration file from the ``cmake.in`` directory, suppose we
are going to use the ``ACML`` library that provides ``BLAS`` and
``LAPACK``.

.. code-block:: bash

    # get the appropriate configuration file

    cp cmake.in/abcdCmake.in.ACML ./abcdCmake.in


To use ``MKL`` instead, copy the file ``abcdCmake.in.MKL``:

.. code-block:: bash

    # get the appropriate configuration file

    cp cmake.in/abcdCmake.in.MKL ./abcdCmake.in

You can use the
`Intel® Math Kernel Library Link Line
Advisor <https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor>`_
to customize the configuration.

Edit the file ``abcdCmake.in`` so that it reflects your configuration (path to libraries, file 
names, path to MPI, etc).


Building the library
--------------------
          
The build process is done using ``cmake``:

.. code-block:: bash

   # create a building directory

   mkdir build

   # run cmake

   cd build
   cmake ..

   # if everything went correctly you can run make

   make

   # the files will be in directory lib/

   ls lib # gives libabcd.a


If cmake does not finish correctly, here are some possible reasons:

* ``mpic++`` is either not installed or there is an issue with ``mpi`` libraries, check also that you 
gave the right path in your ``abcdCmake.in`` file.
* ``Boost`` is either not installed, or the version is too old. Check that ``Boost::MPI`` is installed.
* The path to some libraries or headers is not well defined in ``abcdCmake.in``.

Running ABCD
------------

You can run the solver without having to write a code (as we do in the next section). After building 
the library, a binary is created called ``abcd_run``, it uses a configuration file that you will find in 
the directory ``test/src/config_file.info`` that you need to copy to your build directory.

.. code-block:: bash

   cd build
   cp ../config_file.info .
   
   # to try ABCD on a provided small test matrix, without having to write any code,
   # abcd_run looks by default for the file config_file.info in the current directory

   mpirun -np 16 ./abcd_run

You can also give the executable the path to your configuration file:

.. code-block:: bash

   mpirun -np 16 ./abcd_run /path/to/configuration_file

The configuration file incorporates comments with details about all possible options and how to use them. 
  

Building an example (to call ABCD from C++ or C)
-------------------------------------------------

Once the library is built, you can compile the given examples (either C++ or C):

.. code-block:: bash

   # the C++ example called `example.cpp` and the
   # C example called `example.c` are in the examples directory

   cd examples

   # create a directory where to build your examples

   mkdir build_example
   cd build_example

   # tell cmake where the abcd solver is located
   # the current version supposes that the library was built within
   # the directory ``build`` in a release mode
   # if you get an error while running cmake, check that you gave the
   # absolute path to the abcd solver directory

   cmake .. -DABCD=/absolute/path/to/abcd/
   make

   # if everything went correctly, try to run the C++ example

   mpirun -np 16 ./example

   # or if you want to run the C example:

   mpirun -np 16 ./example_c


Issue tracker
-------------
If you find any bug, didn't understand a step in the documentation, or if you
have a feature request, submit your issue on our
`Issue Tracker <https://bitbucket.org/apo_irit/abcd/issues>`_
by giving:

- reproducible steps
- a source code, or a snippet where you call the solver
- a matrix file if possible.


Release Notes
-------------
* ABCD-1.0
1) Bug fixes:
    a. Patoh imbalance variable is changed to double precision variable.
    b. A new stable uniform partitioning algorithm is implemented.
    c. A new stable algorithm for distributing partitions to MPI processors is implemented. (When the number MPI processes is larger than the number of partitions).
    d. Scaling with MUMPS algorithm is now stable for both new and old versions of MUMPS
2) Improvements:
    a. Now output gives more details.
    b. A post row scaling method is available for getting 2-norm of rows equal 1.
    c. When there is no RHS, the new RHS is created using unscaled input coefficient matrix.
* ABCD-1.1
1) Input/Output:
    a. New parameters:
        + config_file.info:
            - Number of iterations for scaling (manual or predetermined)
        + config_file.info_advanced
            - alpha on the Identity of augmented subsystems: force pivoting to counter numerical issues
            - master_def/slave_def/num_overlap/slave_tol/min_comm_weight
    b. Display:
        + added Block Size/MPI/OpenMP
        + Warning on augmentation Cij without scaling
        + Improved memory display (match MUMPS MB max/avg display)
    c. Partition file example/e05r0500.mtx.part5
    d. Filtering explicit non-zeros of input matrix
2) Improvements:
    a. Overlapping partitions:
        + parameter to specify the number of overlapping rows between contiguous partitions (num_overlap)
        + overlapping lines are at beginning/end of partition
    b. Master-slave scheme
        + Enforce master-slave scheme with new parameter: specify a number of additional slaves (slave_tol)
        + MUMPS analysis:
            - Master with no slave keep their complete analysis
            - The symmetric permutation (METIS/SCOTCH/AMD) is kept between the 2 analysis
        + Explicit distribution of MPI:
            - 1 master/1node + try to fit remaining slaves inside this same node as possible
            - parameters master_def and slave_def allows to choose between old and new implementations
    c. New algorithm to distribute slaves in order to balance communications as well as workload (parameter min_comm_weight to specify imbalance on #rows)
3) Installation system:
    a. dependencies:
        + MUMPS no longer included
        + extraction of Boost for an easier installation
    b. clean CMake/configuration files
    c. documentation
