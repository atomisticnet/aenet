FC      = gfortran -c
#FC      = ifort -c
FCFLAGS = -O2

SHRLD    = gcc
SHRFLAGS = -fPIC
SHARED   = -shared
SHRLIB   = -lm -lgfortran

TAR       = tar
AR        = ar
ARFLAGS   = -crusv

#------------------------------ L-BFGS-B ------------------------------#

LBFGSB = Lbfgsb.3.0

OBJLBFGSB = $(LBFGSB)/blas.o    $(LBFGSB)/lbfgsb.o \
            $(LBFGSB)/linpack.o $(LBFGSB)/timer.o

PICLBFGSB = $(LBFGSB)/blas_pic.o    $(LBFGSB)/lbfgsb_pic.o \
            $(LBFGSB)/linpack_pic.o $(LBFGSB)/timer_pic.o

TRGLBFGSB = liblbfgsb.a
SHRLBFGSB = liblbfgsb.so

#----------------------------------------------------------------------#
#                                Rules                                 #
#----------------------------------------------------------------------#

.SUFFIXES: .o .f $(SUFFIXES)
.PHONY : all static shared clean veryclean

all    : $(LBFGSB) $(TRGLBFGSB) $(SHRLBFGSB)
static : $(LBFGSB) $(TRGLBFGSB)
shared : $(LBFGSB) $(SHRLBFGSB)
lbfgsb : $(LBFGSB) $(TRGLBFGSB) $(SHRLBFGSB)

%.o : %.f Makefile
	$(FC) $(FCFLAGS) $< -o $*.o

%_pic.o : %.f
	$(FC) $(FCFLAGS) $(SHRFLAGS) -o $@ $<

$(LBFGSB) : $(LBFGSB).tar.gz
	$(TAR) xfvz $<

$(TRGLBFGSB) : $(OBJLBFGSB) Makefile
	$(AR) $(ARFLAGS) $@ $(OBJLBFGSB)

$(SHRLBFGSB) : $(PICLBFGSB) Makefile
	$(SHRLD) $(SHARED) $(PICLBFGSB) $(SHRLIB) -o $@

clean :
	if [ -d $(LBFGSB) ]; then rm -rf $(LBFGSB); fi
	for f in $(OBJLBFGSB); do if [ -e $$f ]; then rm -f $$f; fi; done
	for f in $(PICLBFGSB); do if [ -e $$f ]; then rm -f $$f; fi; done
	for f in *~; do rm -f $$f; done

veryclean:
	if [ -d $(LBFGSB) ]; then rm -rf $(LBFGSB); fi
	if [ -e $(TRGLBFGSB) ]; then rm -f $(TRGLBFGSB); fi
	for f in $(OBJLBFGSB); do if [ -e $$f ]; then rm -f $$f; fi; done
	for f in $(PICLBFGSB); do if [ -e $$f ]; then rm -f $$f; fi; done
	for f in *~; do rm -f $$f; done
