##
## Makefile for normaliz
##
include Makefile.configuration

LIBSOURCES = $(wildcard libnormaliz/*.cpp)
LIBHEADERS = $(wildcard libnormaliz/*.h)

SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)

.PHONY : default all
default: normaliz
all: normaliz normaliz1

normaliz: $(SOURCES) $(HEADERS) libnormaliz/libnormaliz.a
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz.cpp libnormaliz/libnormaliz.a $(GMPFLAGS) -o normaliz

normaliz1: $(SOURCES) $(HEADERS) $(LIBHEADERS) $(LIBSOURCES)
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz-impl.cpp $(GMPFLAGS) -o normaliz1

normaliz-pg: $(SOURCES) $(HEADERS) $(LIBHEADERS) $(LIBSOURCES)
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) -pg Normaliz-impl.cpp $(GMPFLAGS) -o normaliz-pg


#always go down the directory and let the make there check what has to be done
.PHONY : libnormaliz/libnormaliz.o
libnormaliz/libnormaliz.o:
	$(MAKE) --directory=libnormaliz libnormaliz.o

.PHONY : libnormaliz/libnormaliz.a
libnormaliz/libnormaliz.a:
	$(MAKE) --directory=libnormaliz libnormaliz.a


.PHONY : clean
clean:
	$(MAKE) --directory=libnormaliz clean
	-rm -f normaliz
	-rm -f normaliz?
	-rm -f normaliz-pg

