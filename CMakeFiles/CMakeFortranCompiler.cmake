SET(CMAKE_Fortran_COMPILER "/usr/bin/gfortran")
SET(CMAKE_Fortran_COMPILER_ARG1 "")
SET(CMAKE_Fortran_COMPILER_ID "GNU")
SET(CMAKE_Fortran_PLATFORM_ID "")

SET(CMAKE_AR "/usr/bin/ar")
SET(CMAKE_RANLIB "/usr/bin/ranlib")
SET(CMAKE_COMPILER_IS_GNUG77 1)
SET(CMAKE_Fortran_COMPILER_LOADED 1)
SET(CMAKE_COMPILER_IS_MINGW )
SET(CMAKE_COMPILER_IS_CYGWIN )
IF(CMAKE_COMPILER_IS_CYGWIN)
  SET(CYGWIN 1)
  SET(UNIX 1)
ENDIF(CMAKE_COMPILER_IS_CYGWIN)

SET(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

SET(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

IF(CMAKE_COMPILER_IS_MINGW)
  SET(MINGW 1)
ENDIF(CMAKE_COMPILER_IS_MINGW)
SET(CMAKE_Fortran_COMPILER_ID_RUN 1)
SET(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;f77;F77;f90;F90;for;For;FOR;f95;F95)
SET(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
SET(CMAKE_Fortran_LINKER_PREFERENCE 20)
IF(UNIX)
  SET(CMAKE_Fortran_OUTPUT_EXTENSION .o)
ELSE(UNIX)
  SET(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
ENDIF(UNIX)

# Save compiler ABI information.
SET(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
SET(CMAKE_Fortran_COMPILER_ABI "")
SET(CMAKE_Fortran_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")

IF(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  SET(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
ENDIF()

IF(CMAKE_Fortran_COMPILER_ABI)
  SET(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
ENDIF(CMAKE_Fortran_COMPILER_ABI)

IF(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  SET(CMAKE_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")
ENDIF()

SET(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;quadmath;m;c")
SET(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/usr/lib/gcc/x86_64-linux-gnu/4.7;/usr/lib/x86_64-linux-gnu;/usr/lib;/lib/x86_64-linux-gnu;/lib")