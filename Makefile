##
## Makefile for normaliz
##
CXX = g++
CXXFLAGS += -Wall -Wno-sign-compare
CXXFLAGS += -O3 -funroll-loops -pipe

## use OpenMP?
ifeq ($(OPENMP),no)
 CXXFLAGS += -Wno-unknown-pragmas -DNO_OPENMP
else
 CXXFLAGS += -fopenmp
endif

NORMFLAGS = -static
N64FLAGS = -Dnorm64 $(NORMFLAGS)
NBIGFLAGS = -Dnormbig $(NORMFLAGS)
GMPFLAGS = -lgmpxx -lgmp

LIBSOURCES = libnormaliz.cpp full_cone.cpp integer.cpp cone_dual_mode.cpp lineare_transformation.cpp list_operations.cpp matrix.cpp mode.cpp output.cpp simplex.cpp sublattice_representation.cpp vector_operations.cpp
LIBHEADERS = $(LIBSOURCES:.cpp=.h)

N64OBJ = obj64/libnormaliz.o obj64/full_cone.o obj64/integer.o obj64/cone_dual_mode.o obj64/lineare_transformation.o obj64/list_operations.o obj64/matrix.o obj64/mode.o obj64/output.o obj64/simplex.o obj64/sublattice_representation.o obj64/vector_operations.o
NBIGOBJ = $(subst obj64,objBig,$(N64OBJ))


all: norm64 normbig

obj64/%.o: %.cpp %.h $(LIBHEADERS)
	@mkdir -p obj64
	$(CXX) $(CXXFLAGS) $(N64FLAGS) -c $< -o $@
norm64: Normaliz.cpp Normaliz.h $(N64OBJ)
	$(CXX) $(CXXFLAGS) $(N64FLAGS) libnormaliz.cpp Normaliz.cpp $(GMPFLAGS) -o norm64

objBig/%.o: %.cpp $(LIBHEADERS)
	@mkdir -p objBig
	$(CXX) $(CXXFLAGS) $(NBIGFLAGS) -c $< -o $@
normbig: Normaliz.cpp Normaliz.h $(NBIGOBJ)
	$(CXX) $(CXXFLAGS) $(NBIGFLAGS) libnormaliz.cpp Normaliz.cpp $(GMPFLAGS) -o normbig


libnormaliz.o: $(LIBHEADERS) $(LIBSOURCES)
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) -c libnormaliz.cpp -o libnormaliz.o 
normaliz: Normaliz.cpp Normaliz.h libnormaliz.o
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz.cpp libnormaliz.o $(GMPFLAGS) -o normaliz


clean:
	-rm -rf obj64 objBig
	-rm -f norm64 normbig
	-rm -f libnormaliz.o normaliz

.PHONY : clean all
