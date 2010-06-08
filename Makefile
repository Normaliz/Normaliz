#
# Makefile for normaliz
#
CXX = g++
CXXFLAGS += -Wall -Wno-sign-compare
CXXFLAGS += -O3 -funroll-loops -pipe
CXXFLAGS += -fopenmp

N64FLAGS = -Dnorm64 -static
NBIGFLAGS = -Dnormbig -static
GMPFLAGS = -lgmpxx -lgmp

SOURCES = full_cone.cpp integer.cpp cone_dual_mode.cpp lineare_transformation.cpp list_operations.cpp matrix.cpp mode.cpp Normaliz.cpp output.cpp simplex.cpp sublattice_representation.cpp vector_operations.cpp
HEADERS = $(SOURCES:.cpp=.h)

N64OBJ = obj64/full_cone.o obj64/integer.o obj64/cone_dual_mode.o obj64/lineare_transformation.o obj64/list_operations.o obj64/matrix.o obj64/mode.o obj64/output.o obj64/simplex.o obj64/sublattice_representation.o obj64/vector_operations.o
NBIGOBJ = $(subst obj64,objBig,$(N64OBJ))


all: norm64 normbig

obj64/%.o: %.cpp %.h $(HEADERS)
	@mkdir -p obj64
	$(CXX) $(CXXFLAGS) $(N64FLAGS) -c $< -o $@
norm64: Normaliz.cpp $(N64OBJ)
	$(CXX) $(CXXFLAGS) $(N64FLAGS) Normaliz.cpp $(N64OBJ) -o norm64

objBig/%.o: %.cpp $(HEADERS)
	@mkdir -p objBig
	$(CXX) $(CXXFLAGS) $(NBIGFLAGS) -c $< -o $@
normbig: Normaliz.cpp $(NBIGOBJ)
	$(CXX) $(CXXFLAGS) $(NBIGFLAGS) Normaliz.cpp $(NBIGOBJ) $(GMPFLAGS) -o normbig

clean:
	-rm -rf obj64 objBig
	-rm -f norm64 normbig

.PHONY : clean all
