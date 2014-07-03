============
Installation
============

Obtaining the source code
-------------------------

The ABCD Solver depends on a few external libraries: ``MUMPS 5.0``, ``Sparselib++ (custom)``, ``PaToH``, ``lapack`` and ``Boost::MPI`` version 1.50 or higher.

* ``MUMPS 5.0``: The latest version of ``MUMPS``, it is distributed
  with our solver in the ``lib/mumps/`` directory. The distributed
  version is compiled in ``x86_64`` compilers, an ``i686`` version can
  be distributed on request.
* ``Sparselib++ (custom)``: a modified version of ``SparseLib++`` to
  suits our needs, is also distributed with our solver in the
  ``lib/sparselib`` directory.
* ``PaToH``: Can be downloaded from
  `Ümit V. Çatalyürek <http://bmi.osu.edu/~umit/software.html>`_.
  webpage. The file ``libpatoh.a`` has to be copied into the ``lib/``
  directory and the header `patho.h` has to be copied into the
  ``include`` directory.

The installation can be done by typing the following commands in your terminal

.. code-block:: bash

    # download the latest stable version
    git clone git@gitlab.enseeiht.fr/mohamed.zenadi/abcd.git

    # we need to download PaToH
    cd abcd/lib

    # download the appropriate version of patoh
    # replace the `xyz` by the build version as described on
    # http://bmi.osu.edu/~umit/software.html
    wget http://bmi.osu.edu/~umit/PaToH/xyz.tar.gz
    # extract
    tar xvzf patoh-xyz.tar.gz
    cp build/xyz/libpatoh.a .
    cp build/xyz/patoh.h ../include/
    rm -rf build patoh-xyz.tar.gz
    cd ..

Now that everything is ready, we can compile the solver. To do so, we need a configuration file from the ``cmake.in`` directory, suppose we are going to use the ``ACML`` library (provides ``blas`` and ``lapack``). 

.. code-block:: bash

    # get the appropriate configuration file
    cp cmake.in/abcdCmake.in.ACML ./abcdCmake.in

Edit that file to suite your configuration. Notice that we link
against ``scalapack`` and ``blacs``, these are required libraries by
``MUMPS``.

To use ``MKL``, copy the file ``abcdCmake.in.MKL``:

.. code-block:: bash

    # get the appropriate configuration file
    cp cmake.in/abcdCmake.in.MKL ./abcdCmake.in

In this example we used the non-threaded version of these libraries,
for the ``ACML`` library use the ``_mp`` suffix, for ``MKL`` use the
`Intel® Math Kernel Library Link Line
Advisor <https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor>`_
to obtain the correct set of libraries.

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
- ``mpic++`` is either not installed or there is an issue with ``mpi`` libraries, check also that you gave the right path in your ``abcdCmake.in`` file.
- ``Boost`` is either not installed, or the version is too old. Check also that ``Boost::MPI`` is installed.
- The path to some libraries is not well defined in ``abcdCmake.in``.
  

Building the example
--------------------

Once the library built, you can compile a given example:

.. code-block:: bash

   # the example.cpp file is in the example directory
   cd example

   mkdir build_example
   cd build_example

   # tell cmake where is the abcd solver
   # the current version supposes that the library was built within
   # the directory ``build`` in a release mode
   # if you get an error while running cmake, check that you gave the
   # absolute path to the abcd solver directory
   cmake .. -DABCD=/absolute/path/to/abcd/
   make

   # if everything went alright 
   mpirun -np 16 ./example
