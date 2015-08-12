# Makefile for genEhrhart

COCOA_ROOT=$(HOME)/CoCoA/CoCoALib-0.99538
include $(COCOA_ROOT)/configuration/autoconf.mk

include ../Makefile.configuration

ALL_C=$(wildcard *.C)
ALL_EXECS=$(ALL_C:.C=)

# Rule for compiling and linking C++ code in *.C files
# LDLIBS contains COCOA_LIB and everythink it needs
nmzIntegrate: nmzIntegrate.C nmzIntInput.C nmzIntPoly.C cyclRatFunct.C genEhrhart.C nmzIntegral.C 
	@echo "Compiling `basename $@`"
	$(COMPILE) -o $@ $< ../libnormaliz/libnormaliz.a $(LDLIBS) #-pthread

.PHONY: library
library:
	@(cd $(COCOA_ROOT); $(MAKE) library)

.PHONY: clean
clean:
	@/bin/rm -f $(ALL_EXECS) *.o

$(ALL_EXECS): $(COCOA_LIB)
