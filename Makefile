##
## Makefile for normaliz
##
include Makefile.configuration

default: normaliz

all: normaliz normaliz1

libnormaliz/libnormaliz.o: $(LIBHEADERS) $(LIBSOURCES)
	$(MAKE) --directory=libnormaliz libnormaliz.o

normaliz: $(SOURCES) $(HEADERS) libnormaliz/libnormaliz.o
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz.cpp libnormaliz/libnormaliz.o $(GMPFLAGS) -o normaliz

normaliz1: $(SOURCES) $(HEADERS) $(LIBHEADERS) $(LIBSOURCES)
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz-impl.cpp $(GMPFLAGS) -o normaliz1


clean:
	$(MAKE) --directory=libnormaliz clean
	-rm -f normaliz

.PHONY : default clean all
