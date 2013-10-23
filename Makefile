##
## Makefile for normaliz
##
include Makefile.configuration

LIBSOURCES = $(wildcard libnormaliz/*.cpp)
LIBHEADERS = $(wildcard libnormaliz/*.h)

SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)

.PHONY : default all linknormaliz
default: linknormaliz
all: lib normaliz normaliz1

linknormaliz: lib
	@$(MAKE) normaliz

normaliz.o: $(SOURCES) $(HEADERS) $(LIBHEADERS)
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) -c Normaliz.cpp -o normaliz.o

normaliz: $(SOURCES) $(HEADERS) normaliz.o libnormaliz/libnormaliz.a
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) normaliz.o libnormaliz/libnormaliz.a $(GMPFLAGS) -o normaliz

normaliz1: $(SOURCES) $(HEADERS) $(LIBHEADERS) $(LIBSOURCES)
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) Normaliz-impl.cpp $(GMPFLAGS) -o normaliz1

normaliz-pg: $(SOURCES) $(HEADERS) $(LIBHEADERS) $(LIBSOURCES)
	$(CXX) $(CXXFLAGS) $(NORMFLAGS) -pg Normaliz-impl.cpp $(GMPFLAGS) -o normaliz-pg


#always go down the directory and let the make there check what has to be done
.PHONY : lib
lib:
	$(MAKE) --directory=libnormaliz libnormaliz.a


.PHONY : clean
clean:
	$(MAKE) --directory=libnormaliz clean
	-rm -f normaliz
	-rm -f normaliz?
	-rm -f normaliz-pg

