##
## Makefile for normaliz
##
CXX = g++
CXXFLAGS += -Wall -Wno-sign-compare -pedantic
CXXFLAGS += -O3 -funroll-loops
#CXXFLAGS += -g #-p

## use OpenMP?
ifeq ($(OPENMP),no)
 CXXFLAGS += -Wno-unknown-pragmas -DNO_OPENMP
else
 CXXFLAGS += -fopenmp
endif

NORMFLAGS = -static
GMPFLAGS = -lgmpxx -lgmp

LIBSOURCES = libnormaliz.cpp cone.cpp full_cone.cpp integer.cpp cone_dual_mode.cpp lineare_transformation.cpp list_operations.cpp matrix.cpp simplex.cpp sublattice_representation.cpp vector_operations.cpp
LIBHEADERS = $(LIBSOURCES:.cpp=.h)

default: libnormaliz.o

libnormaliz.o: $(LIBHEADERS) $(LIBSOURCES) libnormaliz-impl.cpp
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) -c libnormaliz-impl.cpp -o libnormaliz.o 

clean:
	-rm -f libnormaliz.o

.PHONY : default clean all
