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

    # we need to download PaToH
    cd /tmp

    # download the appropriate version of patoh
    # replace the `xyz` by the build version as described in
    # http://bmi.osu.edu/~umit/software.html
    wget http://bmi.osu.edu/~umit/PaToH/xyz.tar.gz
    # extract patoh
    tar xvzf patoh-xyz.tar.gz
    # copy the files to the abcd solver directories lib and include
    cp build/xyz/libpatoh.a /path/to/abcd/lib/
    cp build/xyz/patoh.h /path/to/abcd/include/

    # return to /path/to/abcd/
    cd -

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

   # the files will be in the Release/lib directory
   ls Release/lib # gives libabcd.a


If cmake does not finish correctly, here are some possible reasons:

* ``mpic++`` is either not installed or there is an issue with ``mpi`` libraries, check also that you gave the right path in your ``abcdCmake.in`` file.
* ``Boost`` is either not installed, or the version is too old. Check that ``Boost::MPI`` is installed.
* The path to some libraries is not well defined in ``abcdCmake.in``.
  

Building the example
--------------------

Once the library built, you can compile the given example:

.. code-block:: bash

   # the example.cpp file is in the example directory
   cd example

   mkdir build_example
   cd build_example

   # tell cmake where the abcd solver is located
   # the current version supposes that the library was built within
   # the directory ``build`` in a release mode
   # if you get an error while running cmake, check that you gave the
   # absolute path to the abcd solver directory
   cmake .. -DABCD=/absolute/path/to/abcd/
   make

   # if everything went correctly, try 
   mpirun -np 16 ./example


Issue tracker
-------------

