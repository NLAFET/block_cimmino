============
Installation
============

The ABCD Solver depends on a few libraries: ``MUMPS 5.0``, ``Sparselib++ (custom)``, ``PaToH``, ``lapack`` and ``Boost::MPI``.

The installation can be done by writing the following commands:

.. code-block:: bash

    # download the latest stable version
    git clone https://gitlab.enseeiht.fr/mohamed.zenadi/abcd.git

    # we need to download PaToH
    cd abcd/lib

    # download the appropriate version of patoh
    # replace the ... by the build version as described on
    # http://bmi.osu.edu/~umit/software.html
    wget http://bmi.osu.edu/~umit/PaToH/...tar.gz
    # extract
    tar xvzf patoh-....tar.gz
    cp build/.../libpatoh.a .
    cp build/.../patoh.h ../include/
    rm -rf build patoh-....tar.gz
    cd ..

* ``MUMPS 5.0``: The latest version of ``MUMPS``, it is distributed
  with our solver in the ``lib/mumps/`` directory. The distributed
  version is compiled in ``x86_64`` compilers, a ``i686`` version can
  be distributed on request.
* ``Sparselib++ (custom)``: a modified version of ``SparseLib++`` to
  suits our needs, is also distributed with our solver in the
  ``lib/sparselib`` directory.

- ``PaToH``: Can be downloaded from
  `Ümit V. Çatalyürek <http://bmi.osu.edu/~umit/software.html>`_.
  webpage. The file ``libpatoh.a`` has to be copied into the ``lib/``
  directory and the header `patho.h` has to be copied into the
  ``include`` directory.

Now that everything is ready, we can compile the solver. To do so, we need a configuration file from the ``cmake.in`` directory, suppose we are going to use the ``ACML`` library (provides ``blas`` and ``lapack``). 

.. code-block:: bash

    # get the appropriate configuration file
    cp cmake.in/abcdCmake.in.ACML ./abcdCmake.in

Edit that file to suite your configuration

.. literalinclude:: ../../cmake.in/abcdCmake.in.ACML
                    :language: cmake

Notice that we link against ``scalapack`` and ``blacs``, these are
required libraries by ``MUMPS``. If we want to use ``MKL``, you can
use the file ``cmake.in/abcdCmake.in.MKL``, or just change this part:

.. code-block:: cmake

    set(BLAS_LAPACK_SCALAPACK_DIRS /path/to/mkl)
    set(BLAS_LAPACK_SCALAPACK_LIBS mkl_lapack95_lp64 mkl_blas95_lp64
        mkl_scalapack_lp64 mkl_cdft_core mkl_gf_lp64 mkl_sequential mkl_core
        mkl_blacs_openmpi_lp64 mkl_lapack95_lp64 mkl_blas95_lp64
        )

In this example we used the non-threaded version of these libraries,
for the ``ACML`` library use the ``_mp`` suffix, for ``MKL`` use the
`Intel® Math Kernel Library Link Line
Advisor <https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor>`_
to obtain the correct set of libraries.
