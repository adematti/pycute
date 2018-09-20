####################################################
###          User definable stuff

DEFINEOPTIONS = -D_VERBOSE
#DEFINEOPTIONS += -D_DEBUG
#DEFINEOPTIONS += -D_FLOAT32
#GSL options
GSL_INC = /opt/rhel-6.x86_64/gnu4.6.2/gsl/1.15/include
GSL_LIB = /opt/rhel-6.x86_64/gnu4.6.2/gsl/1.15/lib
CUTE_SRC = ./lib
CUTE_LIB = ./lib

### End of user-definable stuff
####################################################

LGSL = -L$(GSL_LIB) -lgsl -lgslcblas

# DEFINES for the OpenMP version
DEFINEFLAGSCPU = $(DEFINEOPTIONS)

# COMPILER AND OPTIONS
COMPCPU = gcc
COMPGPU = nvcc
OPTCPU = -Wall -O3 -fopenmp $(DEFINEFLAGSCPU)
OPTCPU_GPU = -Wall -O3 $(DEFINEFLAGSGPU)
OPTGPU = -O3 $(DEFINEFLAGSGPU) -arch compute_20 $(OPT_PRECISION) -Xcompiler -Wall

#INCLUDES AND LIBRARIES
INCLUDECOM = -I./lib -I$(GSL_INC)
LIBCPU = $(LGSL) -lm
LIBGPU = $(LGSL) -L$(CUDADIR)/lib64 -lcudart -lpthread -lm

#.so FILES
DEF = $(CUTE_LIB)/define.so
COM = $(CUTE_LIB)/common.so
CORR = $(CUTE_LIB)/correlator.so
BOX = $(CUTE_LIB)/boxes.so
CUTE = $(CUTE_LIB)/cute.so
OFILES = $(DEF) $(COM) $(CORR) $(BOX) $(CUTE)

#.c FILES
CDEF = $(CUTE_LIB)/define.c
CCOM = $(CUTE_LIB)/common.c
CCORR = $(CUTE_LIB)/correlator.c
CBOX = $(CUTE_LIB)/boxes.c
CCUTE = $(CUTE_LIB)/cute.c
IFILES = $(CDEF) $(CCOM) $(CCORR) $(CBOX) $(CCUTE)


#FINAL GOAL
EXE = CUTE

#RULES
default : $(EXE)
#RULE TO MAKE .o's FROM .c's
$(DEF) : $(CDEF) Makefile
	$(COMPCPU) $(OPTCPU) -fPIC -o $@ $(INCLUDECOM) $(LIBCPU) -shared $(IFILES)
$(COM) : $(CCOM) Makefile
	$(COMPCPU) $(OPTCPU) -fPIC -o $@ $(INCLUDECOM) $(LIBCPU) -shared $(IFILES)
$(CORR) : $(CCORR) Makefile
	$(COMPCPU) $(OPTCPU) -fPIC -o $@ $(INCLUDECOM) $(LIBCPU) -shared $(IFILES)
$(BOX) : $(CBOX) Makefile
	$(COMPCPU) $(OPTCPU) -fPIC -o $@ $(INCLUDECOM) $(LIBCPU) -shared $(IFILES)
$(CUTE) : $(CCUTE) Makefile
	$(COMPCPU) $(OPTCPU) -fPIC -o $@ $(INCLUDECOM) $(LIBCPU) -shared $(IFILES)

#RULES TO MAKE THE FINAL EXECUTABLES
$(EXE) : $(OFILES) Makefile

#CLEANING RULES
clean :
	rm -f $(CUTE_LIB)/*.so

cleaner :
	rm -f $(CUTE_LIB)/*.so $(CUTE_LIB)/*~ *~ $(EXE)
