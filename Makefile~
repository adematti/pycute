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

#.c FILES
CDEF = $(CUTE_LIB)/define.c
CCOM = $(CUTE_LIB)/common.c
CCORR = $(CUTE_LIB)/correlator.c
CBOX = $(CUTE_LIB)/boxes.c
CCUTE = $(CUTE_LIB)/cute.c

#.o FILES
ODEF = $(CUTE_LIB)/define.o
OCOM = $(CUTE_LIB)/common.o
OCORR = $(CUTE_LIB)/correlator.o
OBOX = $(CUTE_LIB)/boxes.o
OCUTE = $(CUTE_LIB)/cute.o
OFILES = $(ODEF) $(OCOM) $(OCORR) $(OBOX) $(OCUTE)

#FINAL GOAL
EXE = CUTE

#RULES
default : $(EXE)
#RULE TO MAKE .o's FROM .c's
$(ODEF) : $(CDEF) Makefile
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM) $(LIBCPU)
$(OCOM) : $(CCOM) Makefile
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM) $(LIBCPU)
$(OCORR) : $(CCORR) Makefile
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM) $(LIBCPU)
$(OBOX) : $(CBOX) Makefile
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM) $(LIBCPU)
$(OCUTE) : $(CCUTE) Makefile
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM) $(LIBCPU)

#RULES TO MAKE THE FINAL EXECUTABLES
$(EXE) : $(OFILES) Makefile
	$(COMPCPU) $(OPTCPU) -fPIC -shared $(OFILES) -o cute.so

#CLEANING RULES
clean :
	rm -f $(CUTE_LIB)/*.o
	rm -f $(CUTE_LIB)/*.so
