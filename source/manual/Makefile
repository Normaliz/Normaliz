##
## Makefile configuration for normaliz
##

#######  customize  ######

LIBNORMALIZ_DIR = ../libnormaliz
INCLUDE_DIR = ..
MY_PROG = maxsimplex

##### end customize ######


CXX = g++ 
#CXX = linux32 g++ -m32 -march=i686  #compile it for linux32
CXXFLAGS += -Wall -pedantic -Wno-long-long
CXXFLAGS += -O3 -funroll-loops
CXXFLAGS += -g       ## debugging
#CXXFLAGS += -pg      ## profiling
#CXXFLAGS += -DNDEBUG ## disables asserts


## use OpenMP?
ifeq ($(OPENMP),no)
 CXXFLAGS += -Wno-unknown-pragmas
else
 CXXFLAGS += -fopenmp
# CXXFLAGS += -openmp_report2
endif

## for distributing the executables link static (does not work for mac)
CXXFLAGS += -static 
## for almost static compilation on Mac use
#CXXFLAGS += -static-libgcc
## make it compatible with older Mac versions
#CXXFLAGS += -mmacosx-version-min=10.5

CXXFLAGS += -I$(INCLUDE_DIR)

LINKFLAGS += -L$(LIBNORMALIZ_DIR)
LINKFLAGS += -lnormaliz -lgmpxx -lgmp 
#GMPFLAGS = -I $(GMPDIR) $(GMPDIR).libs/libgmpxx.a $(GMPDIR).libs/libgmp.a #-liomp5 -lpthread
#GMPFLAGS_MIC = -I $(GMPDIR_MIC) $(GMPDIR_MIC).libs/libgmpxx.a $(GMPDIR_MIC).libs/libgmp.a #-liomp5 -lpthread

default: $(MY_PROG) 

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LINKFLAGS)

#%: %.c
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LINKFLAGS)
#
#%.o: %.cpp %.h
#	$(CXX) $(CXXFLAGS) -c -o $@ $<
#
#%.o: %.c %.h
#	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(MY_PROG) *.o
