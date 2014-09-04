============
Installation
============

Obtaining the source code
-------------------------

The ABCD Solver depends on a few external libraries: ``MUMPS``, ``Sparselib++ (custom)``, ``PaToH``, ``lapack`` and ``Boost::MPI`` version 1.50 or higher.

* A patched version of ``MUMPS`` is distributed with our solver in the
  ``lib/mumps/`` directory. Only the headers and a compiled version
  (``Linux x86_64``, other version will be available uppon request) is
  distributed. When ``MUMPS 5.0`` is released, it should be used
  instead.
* ``Sparselib++ (custom)``: a modified version of ``SparseLib++`` to
  suits our needs, is also distributed with our solver in the
  ``lib/sparselib`` directory. The library is compiled same as MUMPS,
  but you still can recompile it by running ``make all`` in
  ``lib/sparselib`` directory.
* ``PaToH``: Can be downloaded from the webpage of `Ümit V. Çatalyürek
  <http://bmi.osu.edu/~umit/software.html>`_ (URL available in the
  following script). The file ``libpatoh.a`` has to be copied into the
  ``lib/`` directory and the header `patho.h` has to be copied into
  the ``include`` directory.
* ``BLAS`` and ``LAPACK`` are both mandatory. We provide
  configurations to build the solver using ``ACML`` and ``MKL``.
* ``BLACS`` and ``ScaLAPACK`` are required by ``MUMPS``, therefore
  they are needed when you link your software with the solver. We
  explicitly require them so that we can build the examples.
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
    # copy libpatoh.a to the lib/ directory
    # copy patoh.h to the include/ directory

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

Edit the file ``abcdCmake.in`` so that it reflects your configuration (path to libraries, file names, path to MPI, etc).


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

* ``mpic++`` is either not installed or there is an issue with ``mpi`` libraries, check also that you gave the right path in your ``abcdCmake.in`` file.
* ``Boost`` is either not installed, or the version is too old. Check that ``Boost::MPI`` is installed.
* The path to some libraries is not well defined in ``abcdCmake.in``.

Running ABCD
------------

You can run the solver without having to write a code (as we do in the next section). After building the library, a binary is created called ``abcd_run``, it uses a configuration file that you will find in the directory ``test/src/config_file.info`` that you need to copy to your build directory.

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
