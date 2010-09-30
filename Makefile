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

LIBSOURCES = $(wildcard libnormaliz/*.cpp)
LIBHEADERS = $(wildcard libnormaliz/*.h)

SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)


default: normaliz

all: normaliz

#libnormaliz/libnormaliz.o: $(LIBHEADERS) $(LIBSOURCES)
#	$(make) $(CXXFLAGS) $(NORMFLAGS) -c libnormaliz.cpp -o libnormaliz.o 
#normaliz: Normaliz.cpp Normaliz.h output.h output.cpp libnormaliz.o
#	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz.cpp libnormaliz.o $(GMPFLAGS) -o normaliz
#	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz.cpp $(GMPFLAGS) -o normaliz


normaliz: $(SOURCES) $(HEADERS) $(LIBHEADERS) $(LIBSOURCES)
#	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz.cpp libnormaliz.o $(GMPFLAGS) -o normaliz
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz.cpp $(GMPFLAGS) -o normaliz


clean:
	-rm -f libnormaliz/libnormaliz.o normaliz

.PHONY : default clean all
