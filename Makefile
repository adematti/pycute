####################################################
###          User definable stuff

EXEC = pycute/cute.so

DEFINEOPTIONS =
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

#INCLUDES AND LIBRARIES
INCLUDECOM = -I./lib -I$(GSL_INC)

#.c FILES
CCOM = $(CUTE_LIB)/common.c
CCORR = $(CUTE_LIB)/correlator.c
CBOX = $(CUTE_LIB)/boxes.c
CCUTE = $(CUTE_LIB)/cute.c

SRC = $(CCOM) $(CCORR) $(CBOX) $(CCUTE)

#.o FILES

OBJ := $(SRC:.c=.o)

#RULES

default: $(EXEC)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(COMPCPU) $(OPTCPU) -fPIC -shared $^ -o $@

%.o: %.c
	$(COMPCPU) $(OPTCPU) -fPIC -o $@ -c $< $(INCLUDECOM) $(LIBCPU)

#CLEANING RULES
clean :
	rm -f $(CUTE_LIB)/*.o
	rm -f $(EXEC)
