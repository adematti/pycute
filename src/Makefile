####################################################
###          User definable stuff
CC ?= gcc
OMPFLAG ?= -fopenmp
LIBDIR ?= ../lib

DEFINEOPTIONS =
#DEFINEOPTIONS += -D_FLOAT32
# COMPILER AND OPTIONS
OPT ?= $(OMPFLAG) $(DEFINEOPTIONS) -pedantic -Wall -Wl,--allow-multiple-definition -Wextra -O3

#.c FILES
SRC = common.c correlator.c boxes.c cute.c
#.o FILES
OBJ := $(SRC:.c=.o)

vpath %.so $(LIBDIR)

cute.so: $(OBJ)
	@mkdir -p $(LIBDIR)
	$(CC) $(OPT) -shared $^ *.o -lm -o $(LIBDIR)/$@

%.o: %.c
	$(CC) $(OPT) -fPIC -o $@ -c $<

all: cute.so

clean:
	rm -rf $(LIBDIR)/cute.so
	rm *.o

.PHONY: all clean