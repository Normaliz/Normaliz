##
## Makefile for libnormaliz
##
include ../Makefile.configuration

LIBSOURCES = $(wildcard *.cpp)
LIBHEADERS = $(wildcard *.h)

default: libnormaliz.o

all: default

libnormaliz.o: $(LIBHEADERS) $(LIBSOURCES) libnormaliz-impl.cpp
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) -c libnormaliz-impl.cpp -o libnormaliz.o 

clean:
	-rm -f libnormaliz.o

.PHONY : default clean all
