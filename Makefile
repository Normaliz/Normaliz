#
# Makefile for normaliz
#
CXX = g++
CXXFLAGS += -O3 -Wall -Wno-sign-compare
CXXFLAGS += -pipe -funroll-loops -fopenmp

N32FLAGS = -Dnorm32 -static
N64FLAGS = -Dnorm64 -static
NBIGFLAGS = -Dnormbig -static
GMPFLAGS = -lgmpxx -lgmp

SOURCES = full_cone.cpp integer.cpp lineare_transformation.cpp list_operations.cpp matrix.cpp mode.cpp Normaliz.cpp output.cpp simplex.cpp vector_operations.cpp
HEADERS = $(SOURCES:.cpp=.h)

N64OBJ = obj64/full_cone.o obj64/integer.o obj64/lineare_transformation.o obj64/list_operations.o obj64/matrix.o obj64/mode.o obj64/output.o obj64/simplex.o obj64/vector_operations.o
N32OBJ = $(subst obj64,obj32,$(N64OBJ))
NBIGOBJ = $(subst obj64,objBig,$(N64OBJ))


all: norm32 norm64 normbig

obj64/%.o: %.cpp %.h $(HEADERS)
	@mkdir -p obj64
	$(CXX) $(CXXFLAGS) $(N64FLAGS) -c $< -o $@
norm64: Normaliz.cpp $(N64OBJ)
	$(CXX) $(CXXFLAGS) $(N64FLAGS) Normaliz.cpp $(N64OBJ) -o norm64

obj32/%.o: %.cpp $(HEADERS)
	@mkdir -p obj32
	$(CXX) $(CXXFLAGS) $(N32FLAGS) -c $< -o $@
norm32: Normaliz.cpp $(N32OBJ)
	$(CXX) $(CXXFLAGS) $(N32FLAGS) Normaliz.cpp $(N32OBJ) -o norm32

objBig/%.o: %.cpp $(HEADERS)
	@mkdir -p objBig
	$(CXX) $(CXXFLAGS) $(NBIGFLAGS) -c $< -o $@
normbig: Normaliz.cpp $(NBIGOBJ)
	$(CXX) $(CXXFLAGS) $(NBIGFLAGS) Normaliz.cpp $(NBIGOBJ) $(GMPFLAGS) -o normbig

clean:
	-rm -rf obj64 obj32 objBig
	-rm -f norm32 norm64 normbig
