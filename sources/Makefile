.SUFFIXES:.out .o .s .c .F .f .f90 .F90 .e .r .y .yr .ye .l .p .sh .csh .h

# use these lines to produce a parallel (mpi) version to parallelize the inversion

FLAGS = -g -lgdal -I/usr/include -c -O \
-fbounds-check \
-Wunused-parameter
#-fcheck-array-temporaries \
-Wall \

CPPFLAGS = -g -lm -lgdal -lstdc++ -c

FF90 = mpif90
FF77 = mpif77 -DD
CC = mpicc
CPP = mpic++


OBJECTS = \
modules.o \
na.o \
Ketcham.o \
read_mtx.o \
CFC.o \
troutq.o \
tripack.o \
locategdal.o \
calc_distance.o


.f90.o:
	$(FF90) $(FLAGS) $*.f90

.F90.o:
	$(FF90) $(FLAGS) $*.F90

.f.o:
	$(FF77) $(FLAGS) $*.f
	
.c.o:
	$(CC) $(FLAGS) $*.c

.cpp.o:
	$(CPP) $(CPPFLAGS) $*.cpp

NCFA:$(OBJECTS)
	$(FF90) $(OBJECTS) -lgdal -lm -lstdc++ -I/usr/include -o ../bin/NCFA
	
all: $(OBJECTS)
	make NFCA
	
clean:
	rm *.o
	rm *.mod
