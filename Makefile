##
## Makefile for normaliz
##
include Makefile.configuration

LIBSOURCES = $(wildcard libnormaliz/*.cpp)
LIBHEADERS = $(wildcard libnormaliz/*.h)

SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)

default: normaliz

all: normaliz normaliz1

libnormaliz/libnormaliz.o: $(LIBHEADERS) $(LIBSOURCES)
	$(MAKE) --directory=libnormaliz libnormaliz.o

libnormaliz/libnormaliz.a: $(LIBHEADERS) $(LIBSOURCES)
	$(MAKE) --directory=libnormaliz libnormaliz.a

normaliz: $(SOURCES) $(HEADERS) libnormaliz/libnormaliz.a
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz.cpp libnormaliz/libnormaliz.a $(GMPFLAGS) -o normaliz

normaliz1: $(SOURCES) $(HEADERS) $(LIBHEADERS) $(LIBSOURCES)
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz-impl.cpp $(GMPFLAGS) -o normaliz1


normaliz-pg: $(SOURCES) $(HEADERS) $(LIBHEADERS) $(LIBSOURCES)
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) -pg Normaliz-impl.cpp $(GMPFLAGS) -o normaliz-pg

clean:
	$(MAKE) --directory=libnormaliz clean
	-rm -f normaliz
	-rm -f normaliz?
	-rm -f normaliz-pg

.PHONY : default clean all
