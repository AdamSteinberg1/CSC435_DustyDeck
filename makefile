# Makefile to build lbstime.a library and driver routines
#
# Andrew J. Pounds, Ph.D.
# Departments of Chemistry and Computer Science
# Mercer University
# Spring 2007
#

F77 = gfortran
CC  = g++
CFLAGS = -O3
FFLAGS = -O0
PROF = -pg
GCOV = -fprofile-arcs -ftest-coverage


TIMINGLIBS =  -L./ -llbstime
CLIBS = -lm

OBJS = cputime.o walltime.o

all: dusty dusty_f90 dusty_c lib

cputime.o : cputime.cc
	$(CC) -O3 -c cputime.cc

walltime.o : walltime.cc
	$(CC) -O3 -c walltime.cc

dusty.o : dusty.f
	$(F77) $(FFLAGS) -c dusty.f

dusty_f90.o : dusty_f90.f90
	$(F77) $(FFLAGS) -c dusty_f90.f90

dusty_c.o : dusty_c.cc
	$(CC) $(CFLAGS) -c dusty_c.cc

dusty : dusty.o lib  $(OBJS)
	$(F77) -o dusty dusty.o  $(TIMINGLIBS) -lstdc++

dusty_f90 : dusty_f90.o lib  $(OBJS)
	$(F77) -o dusty_f90 dusty_f90.o  $(TIMINGLIBS) -lstdc++

dusty_c : dusty_c.o lib  $(OBJS)
	$(CC) -o dusty_c dusty_c.o  $(TIMINGLIBS) -lstdc++

# Default Targets for Cleaning up the Environment
clean :
	rm *.o
	rm *.a

pristine :
	rm *.o
	rm *.a
	touch *.cc *.f
	rm dusty
	rm dusty_f90

ctags :
	ctags  *.cc *.f

# Target for making the library

lib: $(OBJS)
	ar -rc liblbstime.a $(OBJS)
	ranlib liblbstime.a
