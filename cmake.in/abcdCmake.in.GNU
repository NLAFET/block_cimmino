set(BOOST_ROOT /usr/include)

set(BLAS_DIR /usr/lib)
set(LAPACK_DIR /usr/lib)
set(SCALAPACK_DIR /usr/lib)
set(SPARSELIB_DIR $ENV{HOME}/projects/sparselib_1_7)
set(MV_DIR $ENV{HOME}/projects/sparselib_1_7/mv)

set(METIS_DIR /usr/lib)
set(SCOTCH_DIR /usr/lib)
set(MUMPS_DIR /usr/lib)

#----DO NOT CHANGE FROM HERE UNLESS YOU KNOW WHAT YOU ARE DOING----
set(LIB_DIRS    ${BLAS_DIR}
                ${LAPACK_DIR}
                ${SCALAPACK_DIR}
                ${SPARSELIB_DIR}/lib
                ${METIS_DIR}
                ${SCOTCH_DIR}
                ${MUMPS_DIR}
)
set(LIBS
    reflapack   #LAPACK
    refblas     #BLAS
    scalapack   #SCALAPCK + BLACS
    parmetis    
    metis
    esmumps
    scotch
    scotcherr
    dmumps
    mumps_common

    # --- DO NOT REMOVE ---
    sparse
    spblas
    mv
)
