module definitions

  type faulttype

! A fault vertical cross-sectional (2D) geometry is defined by a set of points (x,y)
! Its 3D location is defined by two points in the horizontal plane (x1,y1),(x2,y2)
! Those points limit its lateral extent
! In the (x,y) coordinate system (in which the faul vertical geometry is defined)
! x lies to the right of the (x1,y1),(x2,y2) line and z is vertical upwards
! and zero at the surface

! n is the number of points used to define the fault
! x(n) and y(n) are the coordinates of the points
! (x1,y1) and (x2,y2) are the coordinates of the trace of the fault at the surface
! nstep is the number of time intevals used to define the motion story of the fault
! per interval: timestart,timeend a,d velo
! timestart and timeend are the start and end time of the interval
! velo is the velocity across the fault (along the plane of the fault) for the interval

! By convention as one goes along the fault from point 1 to point n
! it is the block to the right that moves

! The sign of velo determines whether the sign of the x-component of velocity

! from that information we compute xs(n-1),ys(n-1), the direction
! of each of the (n-1) segments connecting the n points

! as well as xn,yn the normal to the (x1,y1)-(x2,y2) line in the z=0 plane

  integer n,nstep
  double precision,dimension(:),pointer::x,y
  double precision x1,y1,x2,y2
  double precision,dimension(:),pointer::timestart,timeend,velo
  double precision,dimension(:),pointer::xs,ys
  double precision xn,yn

  end type faulttype

end module definitions

!---------------------------

program testFault

! this program tests the input files for a Pecube run
! it goes through the time steps and topo given in topo_parameters.txt
! and uses the fault information in fault_parameters.txt to produce
! a series of VTK files that can be visualized in MayaVi

! Per time step, it generates a topo file, a velocity file and a series of fault files

use definitions

implicit none

type (faulttype),dimension(:),pointer::fault,faultp
double precision,dimension(:,:,:),allocatable::x,y,z,vx,vy,vz
integer nfault,i,j,k
integer nx,ny,nz,fltflag,nd
integer nstep,istep,jstep
double precision xmin,xmax,ymin,ymax,zmin,zmax
double precision dt,time,xl,yl,zl,xlon1,xlat1,xlon2,xlat2
double precision timesteps(1000)
character iq*1
real*4 param(1024),range(2,1024)

! read in the number of faults

      open (77,file='input/fault_parameters.txt',status='old')
1     read (77,'(a1)') iq
      if (iq.eq.'$' .or. iq.eq.' ') goto 1
      backspace (77)
      read (77,*) nfault
      close (77)

if (nfault.gt.0) allocate (fault(nfault),faultp(nfault))

print*,'Creating topo files'

! gets the geometry of the box, the resolution of the run
! the time steps from the topo_parameters .txt file

! also generates the topo VTK files

call testTopo (nx,ny,nz,xl,yl,zl,timesteps,nstep,xlon1,xlat1,xlon2,xlat2,fltflag)

! reads in the faults gemetry

if (nfault.eq.0) stop

nd=0
call read_in_fault_parameters (fault,nfault,xlon1,xlat1,xlon2,xlat2,xl,yl,zl,timesteps(nstep+1), &
                               nd,range,param)
print*,'done'
! defines the 3D mesh

xmin=0.d0
xmax=xl
ymin=0.d0
ymax=yl
zmin=0.d0
zmax=zl

allocate (x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz))

do k=1,nz
  do j=1,ny
    do i=1,nx
    x(i,j,k)=xmin+(xmax-xmin)*float(i-1)/float(nx-1)
    y(i,j,k)=ymin+(ymax-ymin)*float(j-1)/float(ny-1)
    z(i,j,k)=zmin+(zmax-zmin)*float(k-1)/float(nz-1)
    enddo
  enddo
enddo

! steps through time backwards from now (t=0) calculating the velocity
! at the time steps defined in topo_parameters.txt
! as well as the geometry of the faults

! doing so it creates VTK files

allocate (vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz))

time=timesteps(nstep+1)

print*,'Creating fault and velocity files'

do istep=nstep,1,-1

! finds the velocity field

do k=1,nz
  do j=1,ny
    do i=1,nx
    call find_velocity (x(i,j,k),y(i,j,k),z(i,j,k), &
                    vx(i,j,k),vy(i,j,k),vz(i,j,k), &
                    fault,nfault,time,1,xmin,xmax,ymin,ymax)
    enddo
  enddo
enddo

! creates the VTK files

call VTK (x,y,z,vx,vy,vz,nx,ny,nz,fault,nfault,istep+1,timesteps(nstep+1)-time)

print*,'Step',istep,'done'

! updates the fault geometry

dt=-(timesteps(istep+1)-timesteps(istep))/100
  do jstep=1,100
  time=time+dt
  if (fltflag.eq.1) call move_fault (fault,faultp,nfault,dt,time)
  enddo
enddo

istep=0

do k=1,nz
  do j=1,ny
    do i=1,nx
    call find_velocity (x(i,j,k),y(i,j,k),z(i,j,k), &
                    vx(i,j,k),vy(i,j,k),vz(i,j,k), &
                    fault,nfault,time,1)
    enddo
  enddo
enddo


call VTK (x,y,z,vx,vy,vz,nx,ny,nz,fault,nfault,istep+1,timesteps(nstep+1)-time)

print*,'Step',istep,'done'

deallocate (x,y)
deallocate (vx,vy,vz)
if (nfault.gt.0) deallocate (fault,faultp)

end program testFault

!----------------

subroutine VTK (x,y,z,vx,vy,vz,nx,ny,nz,fault,nfault,istep,time)

use definitions

implicit none

type (faulttype) fault(nfault)
double precision x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz)
double precision vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
integer nx,ny,nz,nfault
integer iunit,i1,i2,i,j,k,istep,itime
character cs*3,cf*1,ct*3
double precision ymin,ymax,time

write(cs,'(i3)') istep-1
if (istep-1.lt.10) cs(1:2)='00'
if (istep-1.lt.100) cs(1:1)='0'

itime=int(time+1.d-10)
write(ct,'(i3)') itime
if (itime.lt.10) ct(1:2)='00'
if (itime.lt.100) ct(1:1)='0'

iunit=30
open(unit=iunit,file='VTK/velo'//cs//'.vtk')
write(iunit,'(a)')'# vtk DataFile Version 3.0'
write(iunit,'(a)')'velocities'
write(iunit,'(a)')'ASCII'
write(iunit,'(a)')'DATASET UNSTRUCTURED_GRID'
write(iunit,'(a7,i10,a6)')'POINTS ',nx*ny*nz,' float'

do k=1,nz
  do j=1,ny
    do i=1,nx
    write(iunit,'(3f16.11)') x(i,j,k),y(i,j,k),z(i,j,k)
    enddo
  enddo
enddo

write(iunit,'(A6, 2I10)') 'CELLS ',(nx-1)*(ny-1)*(nz-1),9*(nx-1)*(ny-1)*(nz-1)
do k=1,nz-1
  do j=1,ny-1
    do i=1,nx-1
    i1=nx*ny*(k-1)+nx*(j-1)+i-1
    i2=i1+nx*ny
    write(iunit,'(9I10)')8,i1,i1+1,i1+nx,i1+nx+1,i2,i2+1,i2+nx,i2+nx+1
    enddo
  enddo
enddo

write(iunit,'(A11, I10)') 'CELL_TYPES ',(nx-1)*(ny-1)*(nz-1)
do k=1,(nx-1)*(ny-1)*(nz-1)
   write(iunit,'(I2)')11 ! octree  (8 nodes)
end do

write(iunit,'(a11,i10)')'POINT_DATA ',nx*ny*nz

write(iunit,'(a)')'VECTORS velo float'
do k=1,nz
  do j=1,ny
    do i=1,nx
    write(iunit,'(3f18.13)') vx(i,j,k),vy(i,j,k),vz(i,j,k)
    enddo
  enddo
enddo

close (iunit)

ymin=minval(y)
ymax=maxval(y)

do i=1,nfault
cf=char(64+i)
iunit=30
open(unit=iunit,file='VTK/fault'//cf//cs//'.vtk')
write(iunit,'(a)')'# vtk DataFile Version 3.0'
write(iunit,'(a)')'surface'
write(iunit,'(a)')'ASCII'
write(iunit,'(a)')'DATASET UNSTRUCTURED_GRID'
write(iunit,'(a7,i10,a6)')'POINTS ',fault(i)%n*2,' float'
  do k=1,fault(i)%n
  write(iunit,'(3f16.11)') fault(i)%x1+fault(i)%xn*fault(i)%x(k), &
                           fault(i)%y1+fault(i)%yn*fault(i)%x(k),fault(i)%y(k)
  enddo
  do k=1,fault(i)%n
  write(iunit,'(3f16.11)') fault(i)%x2+fault(i)%xn*fault(i)%x(k), &
                           fault(i)%y2+fault(i)%yn*fault(i)%x(k),fault(i)%y(k)
  enddo
write(iunit,'(A6, 2I10)') 'CELLS ',fault(i)%n-1,(1+4)*(fault(i)%n-1)
  do k=1,fault(i)%n-1
  write(iunit,'(9I10)')4,k-1,k,k+fault(i)%n,k-1+fault(i)%n
  enddo
write(iunit,'(A11, I10)') 'CELL_TYPES ',fault(i)%n-1
  do k=1,fault(i)%n-1
  write(iunit,'(I2)')9 ! rectangles
  enddo
close(iunit)
enddo

return

end subroutine VTK
