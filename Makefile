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

all: normaliz normalizl

libnormaliz/libnormaliz.o: $(LIBHEADERS) $(LIBSOURCES)
	$(MAKE) --directory=libnormaliz libnormaliz.o

normaliz: $(SOURCES) $(HEADERS) libnormaliz/libnormaliz.o
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz.cpp libnormaliz/libnormaliz.o $(GMPFLAGS) -o normaliz

#normaliz: $(SOURCES) $(HEADERS) $(LIBHEADERS) $(LIBSOURCES)
#	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz-impl.cpp $(GMPFLAGS) -o normaliz


clean:
	$(MAKE) --directory=libnormaliz clean
	-rm -f normaliz

.PHONY : default clean all
