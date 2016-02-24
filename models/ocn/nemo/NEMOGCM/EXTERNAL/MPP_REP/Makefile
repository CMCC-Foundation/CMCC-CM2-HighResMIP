# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include ifort.inc


.DEFAULT:
	-touch $@
all: mpp_domain_decomposition.x
mpp_domain_decomposition.o: mpp_domain_decomposition.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./mpp_domain_decomposition.f90
SRC = ./mpp_domain_decomposition.f90
OBJ = mpp_domain_decomposition.o
clean: neat
	-rm -f .cppdefs $(OBJ) mpp_domain_decomposition.x
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
mpp_domain_decomposition.x: $(OBJ)
	$(LD) $(OBJ) -o mpp_domain_decomposition.x  $(LDFLAGS)
