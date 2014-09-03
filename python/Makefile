TARGET = python_interface

# location of the Python header files

PYTHON = /usr/include/python2.4

BOOST_INC = /usr/include
BOOST_LIB = /usr/lib/x86_64-linux-gnu/

$(TARGET).so: $(TARGET).o
	g++ -shared -Wl,--export-dynamic \
	$(TARGET).o -L$(BOOST_LIB) -lboost_python -lboost_mpi -lboost_serialization\
	-L${HOME}/libs/scalapack/ -lscalapack
	-L${HOME}/libs/blacs/LIB/ -lblacsF77init -lblacsCinit -lblacs -lblacsF77init -lblacsCinit
	-L/usr/lib/ -lblas -llapack

	-L/usr/lib/openmpi/lib/
set(MPI_LIBRARIES 
    mpi
    mpi_f90
    mpi_f77
    mpi_cxx
	-L/usr/lib/python2.7/config-x86_64-linux-gnu -lpython2.7 \
	-o AbcdPy.so

$(TARGET).o: $(TARGET).cpp
	g++ -I$(PYTHON) -I$(BOOST_INC) -I/usr/lib/openmpi/include/ \
		-c $(TARGET).cpp

clean:
	rm -f *.o *.dll *.a *~ core

