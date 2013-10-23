.SUFFIXES:.out .o .s .c .F .f .f90 .F90 .e .r .y .yr .ye .l .p .sh .csh .h

# use these lines to produce a scalar version using the gfortran compiler
# this removes the mpi commands in NA

FLAGS = -g -c -O3 -ftree-vectorize
FF90 = gfortran
FF77 = gfortran -fd-lines-as-comments
CC = gcc

# use these lines to produce a parallel (mpi) version to parallelize the inversion
# Note that some compilers will not allow to have an "include" command preceded by a "d" as an indication
#   of a commented line; if this is the case for your compiler, the routine na.f will not compile properly
#   You will need to edit the na.f file and remove the "d" in the first column of three lines that
#   start with include 'mpif.h'

#FLAGS = -c -O3 -ftree-vectorize
#FF90 = /usr/local/bin/mpif90
#FF77 = /usr/local/bin/mpif77 -fd-lines-as-code
#CC = mpicc

OBJECTS = \
Mad_He.o \
Mad_Trax.o \
Mad_Trax_Zircon.o \
Pecube.o \
calculate_fault_parameters.o \
calculate_ages.o \
create_pecube_in.o \
define_proc.o \
find_distance.o \
find_dt.o \
find_element.o \
find_neighbours.o \
find_temp.o \
find_triangle.o \
find_velo.o \
find_velocity.o \
four1.o \
interpol3d.o \
interpolate.o \
isostatic_rebound.o \
make_matrix.o \
move_fault.o \
na.o \
read_in_fault_parameters.o \
realft.o \
sinft.o \
solve_iterative.o \
tridag.o \
generate_data.o \
KETCH.o

OBJECTS_TEST = \
testFault.o \
calculate_fault_parameters.o \
find_distance.o \
find_velocity.o \
move_fault.o \
read_in_fault_parameters.o \
testTopo.o

OBJECTS_VTK = \
vtk.o

OBJECTS_NA = \
showNA.o \
create_pecube_in.o \
read_in_fault_parameters.o \
calculate_fault_parameters.o \
find_triangle.o

OBJECTS_CASCADE = \
cascade2Pecube.o

.f90.o:
	$(FF90) $(FLAGS) $*.f90

.F90.o:
	$(FF90) $(FLAGS) $*.F90

.f.o:
	$(FF77) $(FLAGS) $*.f

.c.o:
	$(CC) $(FLAGS) $*.c

Pecube:	$(OBJECTS)
	$(FF90) -O $(OBJECTS) -o ../bin/Pecube 

Test: $(OBJECTS_TEST)
	$(FF90) -O $(OBJECTS_TEST) -o ../bin/Test

Vtk: $(OBJECTS_VTK)
	$(FF90) -O $(OBJECTS_VTK) -o ../bin/Vtk

ShowNA:	$(OBJECTS_NA)
	$(FF90) -O $(OBJECTS_NA) -o ../bin/ShowNA

Cascade2Pecube:	$(OBJECTS_CASCADE)
	$(FF90) -O $(OBJECTS_CASCADE) -o ../bin/Cascade2Pecube

all: $(OBJECTS) $(OBJECTS_TEST)  $(OBJECTS_VTK)
	make Pecube
	make Test
	make Vtk
	make ShowNA
	make Cascade2Pecube

clean:
	rm *.o
