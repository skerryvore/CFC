.SUFFIXES:.out .o .s .c .F .f .f90 .F90 .e .r .y .yr .ye .l .p .sh .csh .h

# use these lines to produce a parallel (mpi) version to parallelize the inversion

FLAGS = -g -lshp -lgdal -I/usr/include/gdal -c -O \
#-fbounds-check \
 -fcheck-array-temporaries \
# -Wall \
# -Wunused-parameter

FF90 = mpif90
FF77 = mpif77 -DD
CC = mpicc

OBJECTS = \
modules.o \
FF3D.o \
na.o \
Ketcham.o \
read_mtx.o \
NCFA.o \
shape.o \
dbf.o \
gdal.o \
troutq.o \
tripack.o


.f90.o:
	$(FF90) $(FLAGS) $*.f90

.F90.o:
	$(FF90) $(FLAGS) $*.F90

.f.o:
	$(FF77) $(FLAGS) $*.f
	
.c.o:
	$(CC) $(FLAGS) $*.c

NCFA:$(OBJECTS)
	$(FF90) $(OBJECTS) -lshp -lgdal -I/usr/include/gdal -o ../bin/NCFA
	
all: $(OBJECTS)
	make NFCA
	
clean:
	rm *.o
	rm *.mod
