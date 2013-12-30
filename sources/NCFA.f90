program NFCA
use iso_c_binding
use precision_kind
use CFC
use fwd
use utilities

implicit none
  include "mpif.h"
  integer,parameter :: LU_input=31
  integer,parameter :: LU_param=32
  integer,parameter :: LU_data=33
  integer,parameter :: LU_output=34
  integer,parameter :: LU_dem=35
  integer,parameter :: LU_flt=36
  integer,parameter :: LU_dummy=37
  character(len = 150 ) :: input_file
  character(len = 150 ) :: data_file
  character(len = 200 ) :: param_file  
  character(len = 150 ) :: output_file
  character(len = 150 ) :: dem_file
  character(len = 200 ) :: filename
  character(len = 200 ) :: tiffname
  character(len = 1024) :: line
  character(len = 20  ) :: dummy
  integer :: UTM_PROJECTION_ZONE(2)
  integer :: ndata,icol
  integer :: i,j,k,io,kk,L,M
  integer :: x1, x2, ny, nx
  logical :: SUPPRESS_UTM_PROJECTION
  real(kind=r4),allocatable,dimension(:,:) :: distances
  integer,allocatable,dimension(:,:) :: dem_topo,flt_topo
  integer :: nx_dem,ny_dem
  real(kind=r4) :: misfit
  real(kind=r8) :: rlon,rlat,rx,ry
  real(kind=r4)  :: range(2,1024)
  type (database),dimension(:),allocatable :: sample
  integer :: nd
  integer :: ndt
  integer :: ierr, r1
  integer :: nodata
  logical :: file_exists
  integer :: x_index,y_index
  integer :: nmaps, ncsamp, ncloc
  real(kind=r4), dimension(:),allocatable :: maps 
  complex(kind=r4),dimension(:,:),allocatable :: CA
  complex(kind=r4),dimension(:),allocatable :: CWK
  real(kind=r4),dimension(:),allocatable :: RWK
  real(kind=r4),dimension(:),allocatable :: IWK
  real(kind=r4) :: xllcorner,yllcorner,cellsize
  real(kind=r4),dimension(:),allocatable :: tspan
  real(kind=r4) :: interp1
  real(kind=r4) :: dt
  real(kind=r4) :: mfitmin
  real(kind=r4) :: ulcorner(6), yulcorner, xulcorner
  real(kind=r4) :: model_opt(50)
  real(kind=r4 ) :: cutwave, samp_freq
  integer :: cutlat, cutlon 
  integer :: mopt
  integer :: nproc
  integer :: iproc
  integer :: CIRCLE_DEF
  integer, parameter :: NCLOSEST=2, RADIUS =1, VORONOI=3 
  logical :: lroot  
  logical :: debug
  logical :: FIRST_PATH, BOREHOLE
  logical :: time_sequence(10)
  real(kind=r4) :: t1,t2
  real(kind=r4) :: MTL, AFTA, FTLD(200), LKH 
  real(kind=r4) :: llTL, llFT 
  real(kind=r4), dimension(:), allocatable :: sorted_distances
  integer, dimension(:), allocatable :: idx
  integer :: sampleid
  CHARACTER(LEN=100), DIMENSION(200), TARGET :: stringArray
  TYPE(C_PTR), DIMENSION(200) :: stringPtrs
  integer :: UTM_zone

  ! List of variables used by TRIPACK
  integer, dimension(:), allocatable :: LIST, LPTR, LEND, NEAR, NEXT
  integer :: IER, LNEW
  real(kind=r8),dimension(:), allocatable :: DIST

  ! List of variables used by TROUTQ
  integer :: NCC, LCC(1), LOUT, NAT, NB, NT
  integer,dimension(:), allocatable :: NNABS, NPTR, NPTR1, NABOR, NBNOS  

  common/invert_type/time_sequence  
  common /NA_MPI/iproc,nproc,lroot
  common /NA_misfit_info/mfitmin,mopt,model_opt,sampleID
  
  ! Variable to be sent to the Ketcham C_routine
  real(kind = c_double),parameter :: alo=16.3
  real(kind = c_double) :: kAFTA,kFTLD(200),kMTL, kold
  real(kind = c_float) :: ktime(4),ktemp(4)

  interface 
    subroutine ketch_main(ntime, ketchtime, ketchtemp, alo, final_age, oldest_age, fmean, fdist) bind(C, name="ketch_main_")
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: ntime
      real(kind=c_float) :: ketchtime(*)
      real(kind=c_float) :: ketchtemp(*)
      real(kind=c_double ) :: alo, final_age, oldest_age, fmean
      real(kind=c_double ) :: fdist(200)
    end subroutine
  end interface

  interface
    subroutine create_shape_points(nSHPType, pszfilename, sampx, sampy, nsamp) bind(C,name="CSHPpoints")
      use iso_c_binding
      implicit none
      integer(kind=c_int),value :: nSHPType
      character(kind=c_char), dimension(*) :: pszfilename
      real(kind = c_double), dimension(*) :: sampx, sampy
      integer(kind=c_int),value :: nsamp 
    end subroutine
  end interface  

  interface
    subroutine createDBF(dbffilename, id, nsamp, sampx, sampy, sampz, sampname &
	&, temp, ERate, Geo, ExRate, Ftage) bind(C, name="createDBF")
      use iso_c_binding
      character(kind=c_char), dimension(*) :: dbffilename
      integer(kind=c_int), value :: nsamp
      integer(kind=c_int), dimension(*) :: id
      real(kind = c_double), dimension(*) ::  sampx, sampy, sampz, temp, ERate
      real(kind = c_double), dimension(*) ::  Geo, ExRate, Ftage
      type(c_ptr), dimension(*) :: sampname
    end subroutine
  end interface
  
  interface
    subroutine create_tiff(tiffname, ulcorner, vals, nx, ny, utm, geoid) bind(c, name="create_tiff")
      use iso_c_binding
      implicit none
      character(kind=c_char), dimension(*) :: tiffname
      real(kind = c_double), dimension(6) :: ulcorner
      integer(kind = c_int), dimension(*) :: vals
      integer(kind = c_int), value :: nx, ny, utm
      character(kind=c_char), dimension(*) :: geoid 
    end subroutine
  end interface
  
  debug=.TRUE.

  call MPI_Init(ierr )      
  call MPI_Comm_size( MPI_COMM_WORLD, nproc, ierr ) 
  call MPI_Comm_rank( MPI_COMM_WORLD, iproc, ierr ) 

  if(iproc.eq.0) call CPU_TIME(t1)
    
  ! open parameters file 
  param_file="input/param_file.txt" 
  call clean_input(LU_dummy,param_file,LU_param,"!")
  
  ! Read input parameters from input file.
  read(LU_param,*) run_name
  read(LU_param,*) data_file
  read(LU_param,*) output_file
  read(LU_param,*) cutwave
  read(LU_param,*) CIRCLE_DEF
  select case(CIRCLE_DEF)
    case(RADIUS)
      read(LU_param,*) circle
    case(NCLOSEST)
      read(LU_param,*) ncloc
  end select
  read(LU_param,*) dem_file ! read name of dem file  
  read (LU_param,'(a1024)') line
  icol=scan(line,':')
  if (icol.ne.0) then
    search_geo=.TRUE.
    nd=2*ttpoints
    read (line(1:icol-1),*) range(1,nd)
    read (line(icol+1:1024),*) range(2,nd)
  else
    nd=2*ttpoints-1
    search_geo=.FALSE.
    read (line,*) geotherm
  endif
  read(LU_param,*) dt
  read(LU_param,*) nmaps
  allocate(maps(nmaps))
  do I = 1, nmaps
    read(LU_param,*) maps(I)
  end do

  circle=circle*1000 ! Convert to m 
  
  ! Check if a filtered topography file exists:
  inquire(FILE="./DEM/"//trim(adjustL(dem_file(1:index(dem_file,".txt",.false.)))//"flt"), exist=file_exists)
  
  
  ! Filter the topography if not already done:
  if(.not.file_exists) then
    write(*,*) "File ",trim(adjustL(dem_file(1:index(dem_file,".txt",.false.)))//".flt"), " does not exist"
    write(*,*) "Performing low-pass filtering using FFT"  
  
    ! Open DEM file 
    open(LU_dem,file="./DEM/"//trim(adjustL(dem_file)), status='unknown')    
    read(LU_dem,'(a5,i20)') dummy, nx_dem
    read(LU_dem,'(a5,i20)') dummy, ny_dem

    ! To use FFT on a 2D array the dimensions (ny, nx) of the array must be
    ! integer powers of 2. In Order to satisfy this condition and to avoid wrap
    ! around effects of filtering we applied zero-padding all around the borders
    ! of the two dimensional array storing the data. 
    
    ! Calculate padding
    ! Calculate nearest superior power of 2 for each dimension
    nx = 2**ceiling(log(real(nx_dem,r8))/log(2._r8))
    ny = 2**ceiling(log(real(ny_dem,r8))/log(2._r8))

    ! Allocate the arrays. This is where the program is going to require a lot
    ! of memory...depending on the size of the DEM.
    allocate(dem_topo(nx,ny))
    allocate(CA(nx,ny))
    allocate(flt_topo(nx,ny))

    ! Zeros the arrays.
    dem_topo = 0._r8
    flt_topo = 0

    ! We need to put the data in the middle part of the 2D array.
    ! Calulate coordinates of top-left corner from where we will load the dem
    ! values.
    x1 = floor((nx - nx_dem) / 2._r8) 
    x2 = floor((ny - ny_dem) / 2._r8)

    ! Read the header of the dem file.
    read(LU_dem,'(a9,f20.5)') dummy,xllcorner
    read(LU_dem,'(a9,f20.5)') dummy,yllcorner  
    read(LU_dem,'(a8,f20.5)') dummy,cellsize
    !read(LU_dem,'(a12,f20.5)') dummy,nodata

    do J=1,ny_dem
      read(LU_dem,*) (dem_topo(I+x1,J+x2),I=1,nx_dem)
    enddo

    close(LU_dem)  

    ! Convert array to complex
    CA = CMPLX(dem_topo,0._r4)
    
    ! Do forward transform (see subroutine file for calling procedure)
    if(debug) write(*,*) "Starting Fast Fourier Calculation"
    allocate(IWK(6*max(nx,ny)+150))
    allocate(RWK(6*max(nx,ny)+150))
    allocate(CWK(ny))
    call FFT3D(CA,nx,ny,nx,ny,1,1,IWK,RWK,CWK)
    if(debug) write(*,*) "Done"

    ! Apply a low-pass filter to keep wavelength greater than cutwave
    ! (meters)
    
    ! Calculate wave number cut-off
    ! The FFT divide the sampling frequency into nx index.
    ! Each index corresponds to a fundamental frequency value. 
    ! For a sampling frequency of 1/cellsize, the value is 1/(cellsize * nx)
    ! The index number for a frequency 1/wavecut is (cellsize * nx) / wavecut. 
    cutlon = ceiling((nx  *  cellsize) / cutwave)
    cutlat = ceiling((ny  *  cellsize) / cutwave)


    do I=cutlon,nx
      do J=cutlat,ny
	  CA(I,J)=0
      enddo  
    enddo 
    
    ! Do reverse fourier transform (see subroutine file for calling procedure)
    if(debug) write(*,*) "Starting Inverse Fast Fourier Calculation"
    call FFT3D(CA,nx,ny,nx,ny,1,-1,IWK,RWK,CWK)
    flt_topo= real(CA)
    if(debug) write(*,*) "Done"


    ! Write a simple ASCII file (Need to be improved to create either a Tiff or
    ! a BMP)

    open(71, file="./DEM/Filtered_topo.asc", status="unknown")

    write(line,*) "ncols",     nx_dem
    write(71,'(a)') trim(adjustL(line))
    write(line,*) "nrows",     ny_dem
    write(71,'(a)') trim(adjustL(line))
    write(line,*) "xllcorner", xllcorner
    write(71,'(a)') trim(adjustL(line))
    write(line,*) "yllcorner", yllcorner
    write(71,'(a)') trim(adjustL(line))
    write(line,*) "cellsize",  cellsize
    write(71,'(a)') trim(adjustL(line))

    do I = 1 + x2, ny_dem + x2
      do J = 1 + x1, nx_dem + x1
       line = " "
       write(line,'(i5)') flt_topo(J,I)
       write(71,'(a1,a)',advance='no') " ",trim(adjustL(line))	
      end do
      write(71,*)
    end do

    close(71)
    ! Create the .prj file
    ! The coordinate system is identical to the input topography file so
    ! we just need to copy the information from the topo.prj
    open(71,file="./DEM/"//trim(adjustL(dem_file(1:index(dem_file,".asc",.false.)))//"prj"),status='old')
    read(71,'(a)') line
    close(71)
    open(71,file="./DEM/Filtered_topo.prj",status="Unknown")
    write(71,'(a)') trim(adjustL(line))
    close(71)

    ! Create a tiff file (requires GDAL)
    ! (Still in developement)
     xulcorner = xllcorner
     yulcorner = yllcorner + ny_dem * cellsize
     ulcorner = (/xulcorner,cellsize, 0., yulcorner, 0., -1 *cellsize/)
     tiffname = "./DEM/filtered_topo.tiff"
     UTM_Zone = 33

     call create_tiff(trim(adjustL(tiffname))//C_NULL_CHAR, real(ulcorner,c_double), &
                      & int(flt_topo(1+x1:nx_dem+x1, 1+x2:ny_dem+x2), c_int),&
                      & int(nx_dem, c_int), int(ny_dem, c_int), int(UTM_Zone, c_int),&
                      & "WGS84"//C_NULL_CHAR)

    ! Now we want to save the filtered topography to a direct access file.
    ! The file will then be used to extract the elevation of the samples.
    inquire(iolength=r1) flt_topo(1,1)
    open(72, access="direct", recl=r1,form='unformatted', status='scratch')

    do I = 1, ny_dem
      do J = 1, nx_dem
       write(72,rec=((I - 1) * nx_dem + J))  flt_topo(J + x1, I + x2)	
      end do
    end do

    inquire(iolength=r1) dem_topo(1,1)
    open(73, access="direct", recl=r1,form='unformatted', status='scratch')
    do I = 1, ny_dem
      do J = 1, nx_dem
       write(73,rec=((I - 1) * nx_dem + J))  dem_topo(J + x1, I + x2)	
      end do
    end do

  endif
  
  ! We don't need the array in memory anymore.
  deallocate(dem_topo)
  deallocate(flt_topo)
  deallocate(CA)
  deallocate(IWK)
  deallocate(RWK)
  deallocate(CWK)


  ! open sample file list
  open(LU_data,file='./data/'//trim(adjustL(data_file)), status='unknown')
  
  ! read number of lines in the data file
  ndata=0
  read(LU_data,*) ! header
  do
    ndata=ndata+1
    read(LU_data,*,iostat=io)
    if (io < 0) then
     ndata=ndata-1
     print*,"End of file: ",ndata," samples found"
     exit
    endif
  enddo
  
  !Check if the folder "results" exists
  file_exists=.FALSE.
  inquire(FILE="./results", exist=file_exists)
  if(.not.file_exists) then
    CALL system("mkdir results")
  endif
  
  !Create RUN folder
  file_exists=.FALSE.
  inquire(FILE="./results/"//trim(adjustL(run_name)), exist=file_exists)
  if(.not.file_exists) then
    CALL system("mkdir ./results/"//trim(adjustL(run_name)))
  endif
  
  ! allocate data array
  allocate(sample(ndata))
  allocate(distances(ndata,ndata))  
  
  ! read data
  rewind(LU_data)
  read(LU_data,*) ! header
  do I=1,ndata
    read(LU_data,*)&
    &sample(I)%s1name,&  
    &sample(I)%lon,&
    &sample(I)%lat,&  
    &sample(I)%strat_old,&
    &sample(I)%strat_young,&
    &sample(I)%measured_elevation,&
    &sample(I)%filepath
    sample(I)%id = I
  enddo

  ! open each MTX file and read data
  do I=1,ndata
    filename="./data/"//trim(adjustL(sample(I)%filepath))
    call read_mtx(filename,&
                  &sample(I)%NCOUNTS,&
                  &sample(I)%NTL,&
                  &sample(I)%zeta,&
                  &sample(I)%rhod,&
                  &sample(I)%nd,&
                  &sample(I)%FTage,&
                  &sample(I)%FTage_err,&
                  &sample(I)%MTL,&
                  &sample(I)%MTL_err,&
                  &sample(I)%stdev,&
                  &sample(I)%stdev_err,&
                  &sample(I)%NS,&
                  &sample(I)%NI,&
                  &sample(I)%TL,&
                  &sample(I)%s2name)
  enddo
  
  ! convert to UTM coordinates
  UTM_PROJECTION_ZONE(1)=33
  UTM_PROJECTION_ZONE(2)=7

  SUPPRESS_UTM_PROJECTION=.FALSE.

  open(16, file="./DEM/Elevation_compare.txt",status='unknown')
  write(16,*) "Measured ", "Extracted ", "Filtered ", "Lat ", "Lon ", "XUTM ", "YUTM "
  
  do I=1,ndata
    rlon=real(sample(I)%lon,r8) ! Conversion to double precision
    rlat=real(sample(I)%lat,r8) ! Conversion to double precision
    call ll2utm (rlon, rlat, rx, ry, UTM_PROJECTION_ZONE, 3)
    sample(I)%x=real(rx,r4) ! Conversion to simple precision
    sample(I)%y=real(ry,r4) ! Conversion to simple precision

    x_index=int((sample(I)%x-xllcorner)/cellsize)
    y_index=int(ny_dem - (sample(I)%y-yllcorner)/cellsize)

    ! Get dem elevation and filtered elevation from direct access file.
    read(72, rec=((y_index - 1) * nx_dem + x_index)) sample(I)%flt_elevation    
    read(73, rec=((y_index - 1) * nx_dem + x_index)) sample(I)%dem_elevation
    
    write(16,*) sample(I)%measured_elevation, sample(I)%dem_elevation, sample(I)%flt_elevation, &
              & sample(I)%lat, sample(I)%lon, sample(I)%x, sample(I)%y
  enddo
 
    
  ! calculate distances. Distances are calculated in a horizontal plane.
  do I=1,ndata
    do J=1,ndata
      distances(I,J)=(sample(I)%x-sample(J)%x)**2+(sample(I)%y-sample(J)%y)**2
      distances(I,J)=sqrt(distances(I,J))     
    enddo
  enddo
  
  ! find the neighbours
  sample%nneighbours = 0

  do I = 1, ndata
    sample(I)%neighbours(:) = 0
  end do

  kk = 0
  
  if(CIRCLE_DEF /= VORONOI) then

  do I = 1,ndata

    select case(CIRCLE_DEF)

    case(RADIUS)
      do J = 1,ndata
	if(distances(I,J).lt.circle) then
	  sample(I)%nneighbours = sample(I)%nneighbours + 1
	  kk = sample(I)%nneighbours
	  sample(I)%neighbours(kk)         = J
	  sample(I)%neighbours_ncounts(kk) = sample(J)%ncounts
	  sample(I)%neighbours_zeta(kk)    = sample(J)%zeta
	  sample(I)%neighbours_rhodos(kk)  = sample(J)%rhod

	  do K=1, sample(J)%ncounts
	    sample(I)%neighbours_NS(kk,K) = sample(J)%NS(K)
	    sample(I)%neighbours_NI(kk,K) = sample(J)%NI(K)       
	  enddo 

	  sample(I)%neighbours_ages(kk)     = sample(J)%FTage
	  sample(I)%neighbours_ages_err(kk) = sample(J)%FTage_err
	  sample(I)%neighbours_ntl(kk)      = sample(J)%ntl

	  do K=1, sample(J)%ntl
	    sample(I)%neighbours_TL(kk,K) = sample(J)%TL(K)
	  enddo

	  sample(I)%neighbours_MTL(kk)     = sample(J)%MTL
	  sample(I)%neighbours_MTL_err(kk) = sample(J)%MTL_err
	  sample(I)%neighbours_offsets(kk) = sample(J)%z-sample(J)%flt_elevation
	  sample(I)%neighbours_offsets(kk) = 0._r4
	endif
      end do

    case(NCLOSEST)

      ! Sort the samples by distance from the closest to the farest
      allocate(sorted_distances(ndata))
      allocate(idx(ndata))
      sorted_distances = distances(I,1:ndata)
      call indexx(sorted_distances, idx)

      ! Select the N closest locations (ncloc) (i.e. samples from boreholes are
      ! considered as a single location)
      K = 1
      ncsamp = 1
      do while(ncsamp.ne.ncloc)
	BOREHOLE = .FALSE.
	if(K.gt.1) BOREHOLE = (sample(idx(K))%lat.eq.sample(idx(K-1))%lat.and.&
	  & sample(idx(K))%lon.eq.sample(idx(K-1))%lon)

	sample(I)%nneighbours = sample(I)%nneighbours + 1

	kk = sample(I)%nneighbours
	sample(I)%neighbours(kk)         = idx(K)
	print*, idx(k)
	sample(I)%neighbours_ncounts(kk) = sample(idx(K))%ncounts
	sample(I)%neighbours_zeta(kk)    = sample(idx(K))%zeta
	sample(I)%neighbours_rhodos(kk)  = sample(idx(K))%rhod

	do L=1, sample(J)%ncounts
	  sample(I)%neighbours_NS(kk,L) = sample(idx(K))%NS(L)
	  sample(I)%neighbours_NI(kk,L) = sample(idx(K))%NI(L)       
	enddo 

	sample(I)%neighbours_ages(kk)     = sample(idx(K))%FTage
	sample(I)%neighbours_ages_err(kk) = sample(idx(K))%FTage_err
	sample(I)%neighbours_ntl(kk)      = sample(idx(K))%ntl

	do L=1, sample(J)%ntl
	  sample(I)%neighbours_TL(kk,L) = sample(idx(K))%TL(L)
	enddo

	sample(I)%neighbours_MTL(kk)     = sample(idx(K))%MTL
	sample(I)%neighbours_MTL_err(kk) = sample(idx(K))%MTL_err
	sample(I)%neighbours_offsets(kk) = sample(idx(K))%z-sample(idx(K))%flt_elevation
	sample(I)%neighbours_offsets(kk) = 0._r4

	if(.not.BOREHOLE) ncsamp = ncsamp + 1
	K = K + 1

	deallocate(sorted_distances)
	deallocate(idx)
      end do

    end select
      
    end do
  end if

  if(CIRCLE_DEF == VORONOI) then 
    ! Initialize arrays, LIST, LPTR, LEND, LNEW, NEAR, NEXT, DIST, IER
    IER = 0
    allocate(LIST(6 * ndata - 12))
    allocate(LPTR(6 * ndata - 12))
    allocate(LEND(ndata))
    allocate(NEAR(ndata))
    allocate(NEXT(ndata))
    allocate(DIST(ndata))

    ! Calculate delaunay triangulation using TRIPACK
    call TRMESH(ndata, sample%x, sample%y, LIST, LPTR, LEND,LNEW, NEAR, NEXT, DIST,IER) 

    if (IER /= 0) print*, "Error during triangulation"

    ! Find neigbours of each sample
    NCC = 0
    allocate(NNABS(ndata))
    allocate(NPTR(ndata))
    allocate(NPTR1(ndata))
    allocate(NABOR(12*ndata)) ! Not sure about the optimum size
    allocate(NBNOS(ndata))

    call TROUTQ(NCC, LCC, ndata, sample%x, sample%y, LIST, LPTR, LEND, LOUT, NNABS, NPTR,&
      &  NPTR1, NABOR, NBNOS, NAT, NB, NT)

    do I=1, ndata

      sample(I)%nneighbours = NNABS(I)
      kk = sample(I)%nneighbours

      sample(I)%neighbours(1:kk)       = NABOR(NPTR(I):NPTR1(I))

      do L=1, kk 

	J = sample(I)%neighbours(L) 
	sample(I)%neighbours_ncounts(L) = sample(J)%ncounts
	sample(I)%neighbours_zeta(L)    = sample(J)%zeta
	sample(I)%neighbours_rhodos(L)  = sample(J)%rhod

	do K=1, sample(J)%ncounts
	  sample(I)%neighbours_NS(L,K) = sample(J)%NS(K)
	  sample(I)%neighbours_NI(L,K) = sample(J)%NI(K)       
	enddo 

	sample(I)%neighbours_ages(L)     = sample(J)%FTage
	sample(I)%neighbours_ages_err(L) = sample(J)%FTage_err
	sample(I)%neighbours_ntl(L)      = sample(J)%ntl

	do K=1, sample(J)%ntl
	  sample(I)%neighbours_TL(L,K) = sample(J)%TL(K)
	enddo

	sample(I)%neighbours_MTL(L)     = sample(J)%MTL
	sample(I)%neighbours_MTL_err(L) = sample(J)%MTL_err
	sample(I)%neighbours_offsets(L) = sample(J)%z-sample(J)%flt_elevation
	sample(I)%neighbours_offsets(L) = 0._r4
      end do

    end do


    deallocate(LIST)
    deallocate(LPTR)
    deallocate(LEND)
    deallocate(NEAR)
    deallocate(NEXT)
    deallocate(DIST)
    deallocate(NNABS)
    deallocate(NPTR)
    deallocate(NPTR1)
    deallocate(NABOR)
    deallocate(NBNOS)

  end if

  ! Create an output showing the relationships between samples (Connection from
  ! centroid to samples circle)
  open(62, file="results/Circles_connections.txt", status="unknown")
  do I = 1, ndata
    do J = 1, sample(I)%nneighbours
      if(I.ne.sample(I)%neighbours(J)) write(62, '(i4,a7,i4)') I, " ----->", sample(I)%neighbours(J)
    end do    
  end do
  close(62)

  ! invert neighbours

  do I=1,ndata
    write(*,'(a36,a20,a10,i4,a4,i3,a1)',advance='no') "Inverting neighbours around sample:",&
     & trim(sample(I)%s1name),"with ID: ", I,"(N=",sample(I)%nneighbours,")"
    
    fwd_ages    = 0._r4
    fwd_ndata   = 0_i4
    fwd_MTL     = 0._r4
    fwd_offsets = 0._r4
    fwd_ncounts = 0_i4
    fwd_zeta    = 0._r4
    fwd_rhodos  = 0._r4
    fwd_ntl     = 0_i4
    fwd_NS      = 0_i4
    fwd_NI      = 0_i4
    fwd_TL      = 0._r4
    
    fwd_ages    = sample(I)%neighbours_ages(1:kk)
    fwd_ndata   = sample(I)%nneighbours
    fwd_MTL     = sample(I)%neighbours_MTL(1:kk)
    fwd_offsets = sample(I)%neighbours_offsets(1:kk)
    fwd_ncounts = sample(I)%neighbours_ncounts(1:kk)
    fwd_zeta    = sample(I)%neighbours_zeta(1:kk)
    fwd_rhodos  = sample(I)%neighbours_rhodos(1:kk)
    fwd_ntl     = sample(I)%neighbours_ntl(1:kk)
 
    do J=1,fwd_ndata

      do k=1,fwd_ncounts(J)
       fwd_NS(J,k) = sample(I)%neighbours_NS(J,K)
       fwd_NI(J,k) = sample(I)%neighbours_NI(J,K)
      enddo

      do k=1,fwd_ntl(J)
        fwd_TL(J,k) = sample(I)%neighbours_TL(J,K)
      enddo
    enddo  
    

    ! Prepare the parameter space for NA
    model_id=I
    sampleid=I
    time_sequence=.FALSE.
    !time_sequence(1:3)=.TRUE. ! The first three parameters are part of a time sequence
    range(1,1:(ttpoints-1))   = 0._r4   ! Time
    range(2,1:(ttpoints-1))   = 500._r4 ! Time
    range(1,4:(2*ttpoints-1)) = 0._r4   ! Temperature
    range(2,4:(2*ttpoints-2)) = 200._r4 ! Temperature
    range(2,(2*ttpoints-1))   = 20._r4  ! Temperature
    
    mfitmin = 0.

    ! ATTENTION: NA calls the forward model which update some variables in the
    ! modules.
        
    call na(nd,range)

    write(*,*) "misfit:", mfitmin !, "Neighbours:",sample(I)%neighbours_ages(1:sample(I)%nneighbours)

    sample(I)%misfit                     = mfitmin
    sample(I)%bestpath(1,1:(ttpoints-1)) = model_opt(1:(ttpoints-1))
    sample(I)%bestpath(1,ttpoints)       = 500._r4
    sample(I)%bestpath(2,:)              = model_opt(4:(2*ttpoints-1))

    if(search_geo.eqv..TRUE.) then
      sample(I)%bestgeotherm = model_opt(nd)
    else
      sample(I)%bestgeotherm = geotherm
    endif
    
  enddo  
  
  
  ! find family and best model
  ! We need to calculate a sample level misfit for each possible thermal history
  ! (i.e. for each fammily that the sample is a potential member of)
  ! Misfit is different for each sample as we only need to compare misfits for
  ! different thermal histories for the SAME sample data.

  do I = 1,ndata
    write(*,*) "looking for family and optimum model for sample: ", I
 
    FIRST_PATH = .TRUE.

    do J = 1,ndata
      do K = 1,sample(J)%nneighbours
	if(sample(J)%neighbours(K) .eq. I) then

	  L = 1
	  do M = ttpoints,1,-1
	    ktime(M) = abs(sample(J)%bestpath(1,L) - sample(J)%bestpath(1,ttpoints))
	    ktemp(M) = sample(J)%bestpath(2,L) + sample(J)%neighbours_offsets(K) * sample(J)%bestgeotherm
	    L = L + 1
	  end do

	  call ketch_main(int(ttpoints,c_int),ktime,ktemp,alo,kAFTA,kOld,kMTL,kFTLD)

	  AFTA = real(kAFTA,r4)
	  FTLD = real(kFTLD,r4)
	  MTL  = real(kMTL, r4)

	  ! Calculate the Log-likelihood
	  call loglike_FT(AFTA, MTL, sample(I)%zeta, sample(I)%rhod, sample(I)%NS,&
	                 & sample(I)%NI, sample(I)%NCOUNTS, llFT)
	  call loglike_TL(200, sample(I)%NTL, sample(I)%TL, FTLD, llTL)

	  LKH = llFT + llTL

	  ! If the log-likelihood is greater than the previous history
	  ! Store the history and record the misfit

	  if (LKH.lt.sample(I)%optimum_LKH.or.FIRST_PATH) then
            FIRST_PATH = .FALSE.	    
	    sample(I)%optimum_LKH = LKH
	    sample(I)%optimum_path(1,:) = sample(J)%bestpath(1,:)
	    sample(I)%optimum_path(2,:) = sample(J)%bestpath(2,:)
	    sample(I)%optimum_geotherm  = sample(J)%bestgeotherm
	    sample(I)%family            = 0
	    sample(I)%family            = sample(J)%neighbours
	    sample(I)%nfamily           = sample(J)%nneighbours
	  end if
	endif
      enddo
    enddo
  enddo

  
  ! Create an output showing the relationships between samples (Connection from
  ! centroid to samples family)
  open(62, file="results/Family_connections.txt", status="unknown")
  do I = 1, ndata
    do J = 1, sample(I)%nfamily
      if(I.ne.sample(I)%family(J)) write(62, '(i4,a7,i4)') I, " ----->", sample(I)%family(J)
    end do    
  end do
  close(62)


  ! write general results to output file
  open(LU_output,file = trim(adjustL(output_file)), status='unknown')
  write(LU_output,*) "Sample Name", "Lowest Misfit", "Time points", "Temp. points" 
  do I=1,ndata  
    write(LU_output,*) sample(I)%s2name,sample(I)%Lwmisfit,sample(I)%bestpath(1,:),&
    &sample(I)%bestpath(2,:)
  enddo

  ! Determine temperature at each time-step (for output purpose)
  do I = 1, ndata
    sample(I)%time_rec(1:nmaps) = maps
    do J = 1, nmaps
      sample(I)%temp_rec(J) = interp1(sample(I)%optimum_path(1,:),sample(I)%optimum_path(2,:), 500 - maps(J))
    end do
  end do

  ! Determine cooling rate at each time-step (for output purpose)
  do I = 1, ndata
    do J = 1, nmaps
    ! Locate points in the time-temperature history
    K = locate(sample(I)%optimum_path(1,:), sample(I)%temp_rec(J))
    sample(I)%Erate_rec(J) = (sample(I)%optimum_path(2,K+1) - sample(I)%optimum_path(2,K)) / &
                           & (sample(I)%optimum_path(1,K+1) - sample(I)%optimum_path(1,K)) 
    end do 
  end do


  ! Write a new lat, lon, temperature file at each output time.
  do J = 1, nmaps 
    write(filename,*) int(maps(J))
    open(71, file="results/Temp"//trim(adjustl(filename))//".txt", status="unknown")
    write(71,*) "X Y Temp E dT/dz dz/dt"
    do I = 1, ndata
      write(71,'(2f15.2, 4f5.2)') sample(I)%x, sample(I)%y, sample(I)%temp_rec(J),sample(I)%Erate_rec(J), &
	& sample(I)%optimum_geotherm, sample(I)%Erate_rec(J) / sample(I)%optimum_geotherm	
    end do
    close(71)
  end do

  ! Create a shapefile at each output time

  do J = 1, nmaps

    write(filename,*) int(maps(J))
    ! For now the code just creates a shapefile with 2 points at a random location
    ! I am still working on the routine that has to be written in C

    call create_shape_points(1_c_int, "Temp"//trim(adjustL(filename))//".shp"//achar(0),&
      &  real(sample(:)%x,c_double), real(sample(:)%y, c_double), ndata)


    ! Create the .prj file associated to the shapefile
    ! The coordinate system is identical to the input topography file so
    ! we just need to copy the information from the topo.prj
    open(71,file="./DEM/"//trim(adjustL(dem_file(1:index(dem_file,".asc",.false.)))//"prj"),status='old')
    read(71,'(a)') line
    close(71)
    open(71,file="Temp"//trim(adjustL(filename))//".prj",status="Unknown")
    write(71,'(a)') trim(adjustL(line))
    close(71)

    ! Create a dbf file
    do I = 1, ndata
      stringArray(I) = sample(I)%s1name//C_NULL_CHAR
      stringPtrs(I) = C_LOC(stringArray(I))
    end do

    call createDBF("Temp"//trim(adjustL(filename))//".dbf"//achar(0),&
      & sample(:)%id, ndata,  real(sample(:)%x, c_double),&
      & real(sample(:)%y, c_double), &
      & real(sample(:)%dem_elevation, c_double), stringPtrs(1:ndata),&
      & real(sample(:)%temp_rec(J), c_double), real(sample(:)%ERate_rec(J), c_double),&
      & real(sample(:)%optimum_geotherm,c_double), &
      & real(sample(:)%ERate_rec(J) / sample(:)%optimum_geotherm, c_double), &
      & real(sample(:)%FTage, c_double))

  end do

  call CPU_TIME(t2)
  write(*,*) 'Computing time: ',(t2 - t1)/60 ,'minutes'


  deallocate(sample)
  deallocate(distances)


  ! close files
  close(LU_input)
  close(LU_param)
  close(LU_output)
  call MPI_FINALIZE(ierr)

end program



subroutine forward(nd,NA_param,misfit)
  use precision_kind
  use fwd
  use qsort_c_module
  use utilities
  use iso_c_binding

  implicit none
  integer :: I,J,n,K,nages,ndistrib,nmtl
  integer,parameter :: NPOINTS=4
  integer :: nd
  integer :: numPDFPts
  integer :: kk
  real(kind=r4) :: TIME(4),TEMP(4),ORDERED_TIME_SEQ(4)
  real(kind=r4) :: ftage(500),ftldistrib(500,200),ftldmean(500)  
  real(kind=r4) :: misfit_ftage,misfit_ftldmean,misfit_ftldistrib,total_misfit_ftldistrib 
  real(kind=r4) :: NA_param(nd)
  real(kind=r4) :: NA_param_type(nd)
  real(kind=r4) :: misfit
  real(kind=r4) :: U238_DECAYCT=1.55125E-10
  real(kind=r4) :: theta
  real(kind=r4) :: rhosrhoi
  real(kind=r4) :: LKH_FTGRAIN(NCOUNTMAX)
  real(kind=r4) :: LKH_FTAGE(NSAMPLEMAX)
  real(kind=r4) :: LKH_TL(NTLMAX)
  real(kind=r4) :: LKH_TLDistrib(NSAMPLEMAX) 
  real(kind=r4) :: LKH_sample(NSAMPLEMAX)
  real(kind=r4) :: total_LKH
  real(kind=r4) :: pdfAxis(200)

  real(kind=c_float) :: ketcham_time(4),ketcham_temp(4)
  real(kind=c_double) :: ketcham_ftage,ketcham_ftldistrib(200),ketcham_ftldmean
  real(kind=c_double) :: oldest_age
  real(kind=c_double),parameter :: alo=16.3

  interface 
    subroutine ketch_main(ntime, ketchtime, ketchtemp, alo, final_age, oldest_age, fmean, fdist) bind(C, name="ketch_main_")
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: ntime
      real(kind=c_float) :: ketchtime(*)
      real(kind=c_float) :: ketchtemp(*)
      real(kind=c_double ) :: alo, final_age, oldest_age, fmean
      real(kind=c_double ) :: fdist(200)
    end subroutine
  end interface

  TIME(1:3)=NA_param(1:3)
  TIME(4)=500._r4
  TEMP(1:4)=NA_param(4:7)
  if(search_geo) geotherm=NA_param(8)   
  ftage=0._r4
  ftldistrib=0._r4
  ftldmean=0._r4

  ! This is an attempt to deal with overlying time constraints
  ! The simple idea is to check the order of the time points picked by NA
  ! and assign a very high misfit value (meaning that the solution has to be rejected)
  ! 1) Order the time points
  ! 2) Compare ordered and non ordered arrays
  ! 3) if non-ordered assign high misfit and goto end of the subroutine ( I don't like goto
  ! but we can do that as a temporary solution
  ! 4) if ordered, keep going

  ! We probably won't have to deal with huge arrays but there's no harm to use an efficient
  ! sorting algorithm. Here I chose quicksort (see the module file)
  ! 1)
  ORDERED_TIME_SEQ=TIME(1:NPOINTS)
  call QsortC(ORDERED_TIME_SEQ)
  do I=1,NPOINTS
    if((ORDERED_TIME_SEQ(I)-TIME(I)) /= 0.0_r4) then
      misfit=1e6
      goto 333
    endif
  enddo

  !========== Calculate fission track ages
  LKH_FTAGE=1000000._4            
  do K=1,fwd_ndata
    ! Calculation is done for each sample in the circle as they have different offsets.
    J=1
    do I=NPOINTS,1,-1
      ketcham_time(I)=abs(TIME(J)-TIME(NPOINTS))
      ketcham_temp(I)=TEMP(J)+fwd_offsets(K)*geotherm
      J=J+1
    enddo

    call ketch_main(NPOINTS,ketcham_time,ketcham_temp,real(alo,8),&
      &ketcham_ftage,oldest_age,ketcham_ftldmean,ketcham_ftldistrib)

    ftage(K)=real(ketcham_ftage,4)
    ftldistrib(K,:)=real(ketcham_ftldistrib*100.,4)
    ftldmean(K)=real(ketcham_ftldmean,4)

    if(fwd_ncounts(K).ne.0) then
      call loglike_FT(ftage(K), ftldmean(K), fwd_zeta(K), fwd_rhodos(K),&
	&  fwd_ns(K,:), fwd_ni(K,:), fwd_ncounts(K), LKH_FTAGE(K))
    endif

    if(fwd_ntl(K).ne.0) then
      call loglike_TL(200, fwd_ntl(K), fwd_TL(K,:), ftldistrib(K,:),LKH_TLDistrib(K))
    endif

    ! Calculate misfit age if no count data (to be implemented)
    ! Calculate misfit mean TL if no TL measurements (to be implemented)
    LKH_sample(K)=LKH_TLDistrib(K)+LKH_FTAGE(K)
  enddo

  !Total likelihood is just the sum of the individual Likelihood.
  total_LKH=sum(LKH_sample(1:fwd_ndata)) 
  !print*,total_LKH
  misfit=abs(total_LKH)        

  333 continue

end subroutine forward


!*************************************************************************
!/*
! * Peter Daly
! * MIT Ocean Acoustics
! * pmd@mit.edu
! * 25-MAY-1998
! * 
! Revisions:
!   Jan. 25, 1999 DHG  Port to Fortran 90
!   Mar. 23, 1999 DHG  To add Lewis Dozier's fix to "rr1" calculation 
! * 
! Description:
! * 
! * These routines convert UTM to Lat/Longitude and vice-versa,
! * using the WGS-84 (GPS standard) or Clarke 1866 Datums.
! * 
! * The formulae for these routines were originally taken from
! * Chapter 10 of "GPS: Theory and Practice," by B. Hofmann-Wellenhof,
! * H. Lictenegger, and J. Collins. (3rd ed) ISBN: 3-211-82591-6,
! * however, several errors were present in the text which
! * made their formulae incorrect.
! *
! * Instead, the formulae for these routines was taken from
! * "Map Projections: A Working Manual," by John P. Snyder
! * (US Geological Survey Professional Paper 1395)
! *
! * Copyright (C) 1998 Massachusetts Institute of Technology
! *               All Rights Reserved
! *
! * RCS ID: $Id: convert_datum.c,v 1.2 1998/06/04 20:50:47 pmd Exp pmd $
! */
!
!*************************************************************************
!
subroutine get_grid_zone (longitude, latitude, grid_zone, lambda0)

  IMPLICIT NONE

  real (kind=8) longitude, latitude
  integer       grid_zone(2)
  real (kind=8) lambda0

  integer  zone_long, zone_lat

  real (kind=8) M_PI
  !!!   parameter (M_PI = 3.141592654)

  !-------------------------------------------------------------------------

  m_pi = ACOS (-1.0)

  !  /* Solve for the grid zone, returns the central meridian */

  zone_long = INT ((longitude + 180.0) / 6.0) + 1
  zone_lat = NINT ((latitude + 80.0) / 8.0)
  grid_zone(1) = zone_long
  grid_zone(2) = zone_lat

  !  /* First, let's take care of the polar regions */

  if ((latitude < -80.0) .OR. (latitude > 84.0)) then
    lambda0 = 0.0 * M_PI / 180.0
    return
  endif

  !  /* Now the special "X" grid */

  if (latitude .GT. 72.0 .AND. &
    longitude .GT. 0.0 .AND. longitude .LT. 42.0) then
  if (longitude .LT. 9.0) then
    lambda0 = 4.5 * M_PI / 180.0
  elseif (longitude .LT. 21.0) then
    lambda0 = 15.0 * M_PI / 180.0
  elseif (longitude .LT. 33.0) then
    lambda0 = 27.0 * M_PI / 180.0
  elseif (longitude .LT. 42.0) then
    lambda0 = 37.5 * M_PI / 180.0
  endif
  return
endif

!  /* Handle the special "V" grid */

if (latitude .GT. 56.0 .AND. latitude .LT. 64.0 .AND. &
  longitude .GT. 0.0 .AND. longitude .LT. 12.0) then
if (longitude .LT. 3.0) then
  lambda0 = 1.5 * M_PI / 180.0
elseif (longitude .LT. 12.0) then
  lambda0 = 7.5 * M_PI / 180.0
endif
return
  endif

  !  /* The remainder of the grids follow the standard rule */

  lambda0 = (FLOAT (zone_long - 1) * 6.0 + (-180.0) + 3.0) * M_PI / 180.0

  return
  end

  !*************************************************************************

  subroutine get_lambda0 (grid_zone, lambda0, ierr)

    IMPLICIT NONE

    integer       grid_zone(2)
    real (kind=8) lambda0
    integer       ierr

    integer zone_long
    integer zone_lat
    real (kind=8) latitude, longitude

    real (kind=8) M_PI
    !!!   parameter (M_PI = 3.141592654)

    !---------------------------------------------------------------------------


    m_pi = ACOS (-1.0)

    !/* Given the grid zone, then set the central meridian, lambda0 */

    !/* Check the grid zone format */

    zone_long = grid_zone(1)
    zone_lat = grid_zone(2)
    if ((zone_long .LT. 1) .OR. (zone_long .GT. 61)) then
      write (*,*) 'Invalid grid zone format: ', zone_long, zone_lat
      ierr = -1
      return 
    endif

    longitude = (FLOAT (zone_long - 1) * 6.0) - 180.0
    latitude = (FLOAT (zone_lat) * 8.0) - 80.0

    !/* Take care of special cases */

    if ((latitude .LT. -80.0) .OR. (latitude .GT. 84.0)) then
      lambda0 = 0.0
      ierr = 0
      return 
    endif

    if (latitude .GT. 56.0 .AND. latitude .LT. 64.0 .AND. &
      longitude .GT. 0.0 .AND. longitude .LT. 12.0) then
    if (longitude .LT. 3.0) then
      lambda0 = 1.5 * M_PI / 180.0
    elseif (longitude .LT. 12) then
      lambda0 = 7.5 * M_PI / 180.0
    endif
    ierr = 0
    return
  endif

  if (latitude .GT. 72.0 .AND. &
    longitude .GT. 0.0 .AND. longitude < 42.0) then
  if (longitude .LT. 9.0) then
    lambda0 = 4.5 * M_PI / 180.0
  elseif (longitude .LT. 21.0) then
    lambda0 = 15.0 * M_PI / 180.0
  elseif (longitude .LT. 33.0) then
    lambda0 = 27.0 * M_PI / 180.0
  elseif (longitude .LT. 42.0) then
    lambda0 = 37.5 * M_PI / 180.0
  endif
  ierr = 0
  return
endif

!/* Now handle standard cases */

lambda0 = (FLOAT (zone_long - 1) * 6.0 + (-180.0) + 3.0) * M_PI / 180.0

!/* All done */

ierr = 0
return
end

!*************************************************************************
subroutine ll2utm (longitude, latitude, utm_x, utm_y, grid_zone, datum)

  IMPLICIT NONE

  real (kind=8) latitude, longitude
  real (kind=8) utm_x, utm_y
  integer       grid_zone(2)
  integer       datum

  real (kind=8)  a, b, f, e, e2, e4, e6
  real (kind=8)  phi, lambda, lambda0, phi0, k0
  real (kind=8)  t, rho, m, x, y, k, mm, mm0
  real (kind=8)  aa, aa2, aa3, aa4, aa5, aa6
  real (kind=8)  ep2, nn, tt, cc

  real (kind=8) M_PI
  !!!   parameter (M_PI = 3.141592654)

  integer CLARKE_1866_DATUM
  parameter (CLARKE_1866_DATUM = 1)
  integer GRS_80_DATUM
  parameter (GRS_80_DATUM = 2)
  integer WGS_84_DATUM
  parameter (WGS_84_DATUM = 3)

  !---------------------------------------------------------------------------


  m_pi = ACOS (-1.0)

  !/* Converts lat/long to UTM, using the specified datum */

  if (datum == CLARKE_1866_DATUM) then      ! CLARKE_1866_DATUM:
    a = 6378206.4
    b = 6356583.8
  elseif (datum == GRS_80_DATUM) then      ! GRS_80_DATUM:
    a = 6378137
    b = 6356752.3
  elseif (datum == WGS_84_DATUM) then      ! WGS_84_DATUM:
    a = 6378137.0           !/* semimajor axis of ellipsoid (meters) */
    b = 6356752.31425       !/* semiminor axis of ellipsoid (meters) */
  else
    write (*,*) 'Unknown datum: ', datum
    return
  endif

  !/* Calculate flatness and eccentricity */

  f = 1 - (b / a)
  e2 = 2 * f - f * f
  e = sqrt (e2)
  e4 = e2 * e2
  e6 = e4 * e2

  !/* Convert latitude/longitude to radians */

  phi = latitude * M_PI / 180.0
  lambda = longitude * M_PI / 180.0

  !/* Figure out the UTM zone, as well as lambda0 */

  call get_grid_zone (longitude, latitude, grid_zone, lambda0)

  phi0 = 0.0

  !/* See if this will use UTM or UPS */

  if (latitude .GT. 84.0) then

    !/* use Universal Polar Stereographic Projection (north polar aspect) */

    k0 = 0.994

    t = sqrt ( ((1 - sin (phi)) / (1 + sin (phi))) * &
      (((1 + e * sin (phi)) / (1 - e * sin (phi))) ** e) )
    rho = 2.0 * a * k0 * t / sqrt ( ((1.0 + e) ** (1.0 + e)) * ((1.0 - e) ** (1.0 - e)) )
    !!! Not needed (dhg) m = cos (phi) / sqrt (1.0 - e2 * sin (phi) * sin (phi))

    x = rho * sin (lambda - lambda0)
    y = -rho * cos (lambda - lambda0)
    !!! Not needed (dhg) k = rho * a * m

    !/* Apply false easting/northing */

    x = x + 2000000.0
    y = y + 2000000.0

  elseif (latitude .LT. -80.0) then

    !/* use Universal Polar Stereographic Projection (south polar aspect) */

    phi = -phi
    lambda = -lambda
    lambda0 = -lambda0

    k0 = 0.994

    t = sqrt (((1.0 - sin (phi)) / (1.0 + sin (phi))) * &
      ( ( (1.0 + e * sin (phi)) / (1.0 - e * sin (phi)) ** e) ) )
    rho = 2.0 * a * k0 * t / sqrt ( ((1+e) ** (1+e)) * ((1-e) ** (1-e)) )
    !!! Not needed (dhg) m = cos (phi) / sqrt (1.0 - e2 * sin (phi) * sin (phi))

    x = rho * sin (lambda - lambda0)
    y = -rho * cos (lambda - lambda0)
    !!! Not needed (dhg) k = rho * a * m

    x = -x
    y = -y

    !/* Apply false easting/northing */

    x = x + 2000000.0
    y = y + 2000000.0

  else

    !/* Use UTM */

    !/* set scale on central median (0.9996 for UTM) */

    k0 = 0.9996

    mm = a * ((1.0-e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0) * phi - &
      (3.0*e2/8.0 + 3.0*e4/32.0 + 45.0*e6/1024.0) * sin (2.0*phi) + &
      (15.0*e4/256.0 + 45.0*e6/1024.0) * sin (4.0*phi) - &
      (35.0*e6/3072.0) * sin (6.0*phi))

    mm0 = a * ((1.0-e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0) * phi0 - &
      (3.0*e2/8.0 + 3.0*e4/32.0 + 45.0*e6/1024.0) * sin (2.0*phi0) + &
      (15.0*e4/256.0 + 45.0*e6/1024.0) * sin (4.0*phi0) - &
      (35.0*e6/3072.0) * sin (6.0*phi0))

    aa = (lambda - lambda0) * cos(phi)
    aa2 = aa * aa
    aa3 = aa2 * aa
    aa4 = aa2 * aa2
    aa5 = aa4 * aa
    aa6 = aa3 * aa3

    ep2 = e2 / (1.0 - e2)
    nn = a / sqrt (1.0 - e2 * sin (phi) * sin (phi))
    tt = tan (phi) * tan (phi)
    cc = ep2 * cos (phi) * cos (phi)

    !!! Not needed (dhg) k = k0 * (1 + (1+cc)*aa2/2 + (5-4*tt+42*cc+13*cc*cc-28*ep2) * aa4 / 24.0 + &
    !!! Not needed (dhg) 	     (61-148*tt+16*tt*tt) * aa6 / 720.0)
    x = k0 * nn * (aa + (1-tt+cc) * aa3 / 6 + &
      (5-18*tt+tt*tt+72*cc-58*ep2) * aa5 / 120.0)
    y = k0 * (mm - mm0 + nn * tan (phi) * &
      (aa2 / 2 + (5-tt+9*cc+4*cc*cc) * aa4 / 24.0 + &
      (61 - 58*tt + tt*tt + 600*cc - 330*ep2) * aa6 / 720))

    !/* Apply false easting and northing */

    x = x + 500000.0
    if (y .LT. 0.0) then
      y = y + 10000000.0
    endif
  endif

  !/* Set entries in UTM structure */

  utm_x = x
  utm_y = y

  !/* done */

  return
  end

  !*************************************************************************
  subroutine utm2ll (utm_x, utm_y, longitude, latitude, grid_zone, datum)

    IMPLICIT NONE

    real (kind=8) utm_x, utm_y
    real (kind=8) latitude, longitude
    integer       grid_zone(2)
    integer       datum

    integer ierr
    real (kind=8)  a, b, f, e, e2, e4, e6, e8
    real (kind=8)  lambda0, x, y, k0, rho, t, chi, phi, phi1, phit
    real (kind=8)  lambda, phi0, e1, e12, e13, e14
    real (kind=8)  mm, mm0, mu, ep2, cc1, tt1, nn1, rr1
    real (kind=8)  dd, dd2, dd3, dd4, dd5, dd6

    real (kind=8) M_PI
    !!!   parameter (M_PI = 3.141592654)
    real (kind=8) LOWER_EPS_LIMIT
    parameter (LOWER_EPS_LIMIT = 1.0e-14)
    real (kind=8) M_PI_2

    integer CLARKE_1866_DATUM
    parameter (CLARKE_1866_DATUM = 1)
    integer GRS_80_DATUM
    parameter (GRS_80_DATUM = 2)
    integer WGS_84_DATUM
    parameter (WGS_84_DATUM = 3)

    !---------------------------------------------------------------------------


    m_pi = ACOS (-1.0)

    M_PI_2 = M_PI * 2.0

    !/* Converts UTM to lat/long, using the specified datum */

    if (datum == CLARKE_1866_DATUM) then      ! CLARKE_1866_DATUM:
      a = 6378206.4
      b = 6356583.8
    elseif (datum == GRS_80_DATUM) then       ! GRS_80_DATUM:
      a = 6378137
      b = 6356752.3
    elseif (datum == WGS_84_DATUM) then       ! WGS_84_DATUM:
      a = 6378137.0             !/* semimajor axis of ellipsoid (meters) */
      b = 6356752.31425         !/* semiminor axis of ellipsoid (meters) */
    else
      write (*,*) 'Unknown datum: ', datum
      return
    endif

    !/* Calculate flatness and eccentricity */

    f = 1.0 - (b / a)
    e2 = (2.0 * f) - (f * f)
    e = sqrt (e2)
    e4 = e2 * e2
    e6 = e4 * e2
    e8 = e4 * e4

    !/* Given the UTM grid zone, generate a baseline lambda0 */

    call get_lambda0 (grid_zone, lambda0, ierr)
    if (ierr .NE. 0) then
      write (*,*) 'Unable to translate UTM to LL'
      return
    endif

    latitude = (FLOAT (grid_zone(2)) * 8.0) - 80.0

    !/* Take care of the polar regions first. */

    if (latitude .GT. 84.0) then !/* north polar aspect */

      !/* Subtract the false easting/northing */

      x = utm_x - 2000000.0
      y = utm_y - 2000000.0

      !/* Solve for inverse equations */

      k0 = 0.994
      rho = sqrt (x*x + y*y)
      t = rho * sqrt ( ((1+e) ** (1+e)) * ((1-e) ** (1-e)) ) / (2*a*k0)

      !/* Solve for latitude and longitude */

      chi = M_PI_2 - 2 * atan (t)
      phit = chi + (e2/2 + 5*e4/24 + e6/12 + 13*e8/360) * sin(2*chi) + &
	(7*e4/48 + 29*e6/240 + 811*e8/11520) * sin(4*chi) + &
	(7*e6/120 + 81*e8/1120) * sin(6*chi) + &
	(4279*e8/161280) * sin(8*chi)

      do while (ABS (phi-phit) .GT. LOWER_EPS_LIMIT)
	phi = phit
	phit = M_PI_2 - 2 * atan ( t * (((1 - e * sin (phi)) / (1 + e * sin (phi))) ** (e / 2)) )
      enddo

      lambda = lambda0 + atan2 (x, -y)

    elseif (latitude .LT. -80.0) then !/* south polar aspect */

      !/* Subtract the false easting/northing */

      x = -(utm_x - 2000000)
      y = -(utm_y - 2000000)

      !/* Solve for inverse equations */

      k0 = 0.994
      rho = sqrt (x*x + y*y)
      t = rho * sqrt ( ((1+e) ** (1+e)) * ((1-e) ** (1-e)) ) / (2*a*k0)

      !/* Solve for latitude and longitude */

      chi = M_PI_2 - 2 * atan (t)
      phit = chi + (e2/2 + 5*e4/24 + e6/12 + 13*e8/360) * sin (2*chi) + &
	(7*e4/48 + 29*e6/240 + 811*e8/11520) * sin (4*chi) + &
	(7*e6/120 + 81*e8/1120) * sin (6*chi) + &
	(4279*e8/161280) * sin (8*chi)

      do while (ABS (phi-phit) .GT. LOWER_EPS_LIMIT)
	phi = phit;
	phit = M_PI_2 - 2 * atan (t * ( ((1-e*sin(phi)) / (1+e*sin(phi)) ) ** (e/2)))
      enddo

      phi = -phi
      lambda = -(-lambda0 + atan2 (x , -y))

    else

      !/* Now take care of the UTM locations */

      k0 = 0.9996

      !/* Remove false eastings/northings */

      x = utm_x - 500000.0
      y = utm_y

      if (latitude .LT. 0.0) then  !/* southern hemisphere */
	y = y - 10000000.0
      endif

      !/* Calculate the footpoint latitude */

      phi0 = 0.0
      e1 = (1.0 - sqrt (1.0-e2)) / (1.0 + sqrt (1.0-e2))
      e12 = e1 * e1
      e13 = e1 * e12
      e14 = e12 * e12

      mm0 = a * ((1.0-e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0) * phi0 - &
	(3.0*e2/8.0 + 3.0*e4/32.0 + 45.0*e6/1024.0) * sin (2.0*phi0) + &
	(15.0*e4/256.0 + 45.0*e6/1024.0) * sin (4.0*phi0) - &
	(35.0*e6/3072.0) * sin (6.0*phi0))
      mm = mm0 + y/k0;
      mu = mm / (a * (1.0-e2/4.0-3.0*e4/64.0-5.0*e6/256.0))

      phi1 = mu + (3.0*e1/2.0 - 27.0*e13/32.0) * sin (2.0*mu) + &
	(21.0*e12/16.0 - 55.0*e14/32.0) * sin (4.0*mu) + &
	(151.0*e13/96.0) * sin (6.0*mu) + &
	(1097.0*e14/512.0) * sin (8.0*mu)

      !/* Now calculate lambda and phi */

      ep2 = e2 / (1.0 - e2)
      cc1 = ep2 * cos (phi1) * cos (phi1)
      tt1 = tan (phi1) * tan (phi1)
      nn1 = a / sqrt (1.0 - e2 * sin (phi1) * sin (phi1))
      !!!DHG Old Code rr1 = a * (1.0 - e2) / ((1.0 - e2 * sin (phi) * sin (phi)) ** 1.5)
      !!!DHG L.Dozier's fix is next
      rr1 = a * (1.0 - e2) / ((1.0 - e2 * sin (phi1) * sin (phi1)) ** 1.5)
      dd = x / (nn1 * k0)

      dd2 = dd * dd
      dd3 = dd * dd2
      dd4 = dd2 * dd2
      dd5 = dd3 * dd2
      dd6 = dd4 * dd2

      phi = phi1 - (nn1 * tan (phi1) / rr1) * &
	(dd2/2.0 - (5.0+3.0*tt1+10.0*cc1-4.0*cc1*cc1-9.0*ep2) * dd4 / 24.0 + &
	(61.0+90.0*tt1+298.0*cc1+45.0*tt1*tt1-252.0*ep2-3.0*cc1*cc1) * dd6 / 720.0)
      lambda = lambda0 + &
	(dd - (1.0+2.0*tt1+cc1) * dd3 / 6.0 + &
	(5.0-2.0*cc1+28.0*tt1-3.0*cc1*cc1+8.0*ep2+24.0*tt1*tt1) * dd5 / 120.0) / cos (phi1)
    endif

    !/* Convert phi/lambda to degrees */

    latitude = phi * 180.0 / M_PI
    longitude = lambda * 180.0 / M_PI

    !/* All done */

    return
    end


    subroutine clean_input(U1,name,U2,char)
      implicit none

      character*200,intent(in) :: name
      integer,intent(in) :: U1
      integer,intent(in) :: U2
      character,intent(in) :: char*1
      character :: line*1024

      integer :: i,k,j
      character*9,parameter :: FMT1='(a1024)'   !Length of 'line'

      open (U1,file=trim(adjustL(name)),status='old')
      open (U2,status='scratch')
      1 read (U1,FMT1,end=2) line
      if (line(1:1).ne.char.and. line(1:1).ne.' ') then
	if (scan(line,char).ne.0) then
	  do i=scan(line,char),1024
	    line(i:i)=' '
	  enddo
	endif
	k=1
	do j=1,1024
	  if (line(j:j).eq.' '.or.line(j:j).eq.',') then
	    if (j.ne.k) write (U2,'(a)') line(k:j-1)
	    k=j+1
	  endif
	enddo
      endif
      goto 1
      2 close (U1)
      rewind (U2)

    end subroutine

    function interp1(x,Y,xi)

      use precision_kind
      implicit none
      integer :: I
      real(kind=r4) :: xi
      real(kind=r4) :: interp1
      !Description

      ! yi = interp1(x,Y,xi) interpolates to find yi, the values of the underlying function Y
      ! at the points in the vector or array xi. x must be a vector. Y can be a scalar, a vector,
      ! or an array of any dimension, subject to the following conditions:
      ! If Y is a vector, it must have the same length as x. A scalar value for Y is expanded
      ! to have the same length as x. xi can be a scalar, a vector, or a multidimensional array,
      ! and yi has the same size as xi.

      real(kind=r4) :: x(*)
      real(kind=r4) :: Y(*)

      I=2
      do while(xi.gt.x(I))
	I=I+1
      enddo

      interp1=(Y(I-1)*(x(I)-xi)+Y(I)*(xi-x(I-1)))/(x(I)-x(I-1))

      return
    end function

    !---------------------------------------------------------------------------------------------------

    subroutine user_init (nd,range,scales)

      use precision_kind
      implicit none

      real(kind=r4),intent(in) :: range(2,*)
      real(kind=r4),intent(out) :: scales(*)
      integer,intent(in) :: nd

      scales(1:nd)=-1.

      return
    end subroutine

    !---------------------------------------------------------------------------------------------------
    subroutine writemodels (nd,ntot,models,misfit,ns1,ns2,itmax, &
	nh_max,nh,header)

      use precision_kind
      use fwd

      implicit none

      integer :: nd
      integer :: ntot
      integer :: ns1
      integer :: ns2
      integer :: itmax
      integer :: nh_max
      integer :: nh
      integer :: i
      integer :: k

      real(kind=r4) :: models (nd,*)
      real(kind=r4) :: misfit (ntot)
      character*(*) header
      character*5 :: tempo

      write(tempo,'(i5)') model_id
      open (87,file='./results/'//trim(adjustL(run_name))//'/NA_results_sample'//trim(adjustL(tempo))//'.txt',status='unknown')
      do i=1,ntot
	write (87,*) misfit(i),(models(k,i),k=1,nd)
      enddo
      close (87)

      return
    end subroutine


    !-----------------------------------------------------------------------
    !
    ! subroutine loglike_FT
    !
    ! 
    ! For fission track data, a natural choice of data fit is the log 
    ! likelihood function given by Gallagher (1995). This is defined in terms of the
    ! observed spontaneous and induced track counts, Nsj and Nij for each crystal j of
    ! a total of Nc.   
    !
    subroutine loglike_FT(PredAFTA, PredMTL, zeta, rdos, ns, ni, nc, lkh)
      use precision_kind
      implicit none

      real(kind = r4), intent(in) :: PredAFTA, PredMTL, zeta, rdos
      integer, intent(in) :: nc, ns(nc), ni(nc) 
      real(kind = r4) :: lkh

      integer :: j
      real(kind = r4), parameter :: U238_DECAYCT = 1.55125E-10, alo = 16.3
      real(kind = r4) :: rsri, theta  


      ! 1) First, calculate rhos/rhoi from predicted age
      rsri = (exp(U238_DECAYCT * PredAFTA * 1D6) - 1)
      rsri = rsri * (2 * (PredMTL / alo)) / (U238_DECAYCT * zeta * rdos)

      ! 2) Calculate theta from rhos/rhoi
      theta = rsri / (1 + rsri)

      ! 3) Calculate log-likelihood using count data
      lkh = 0._r4
      do j = 1, nc
	lkh = lkh + ns(j) * log(theta) + ni(j) * log(1 - theta)	
      end do

    end subroutine loglike_FT  

    !-----------------------------------------------------------------------
    !
    ! subroutine loglike_TL
    subroutine loglike_TL(nx, ntl, tl, ftld, lkh)
      use precision_kind
      implicit none

      integer, intent(in) :: nx,ntl 
      real(kind = r4),intent(in) :: tl(ntl), ftld(nx)
      real(kind = r4),intent(out) :: lkh 

      integer :: i, j
      real(kind = r4) :: pdfx(nx) 

      do i = 1, nx
	pdfx(i) = (i-1) * 0.1_r4	
      end do

      lkh = 0._r4

      do i = 1, ntl
	j = locate(pdfx, tl(i))
	lkh = lkh + log(ftld(j) / 100.0)
      end do


    contains


      function locate(xx,x)
	use precision_kind
	implicit none

	real(kind=r4),dimension(:),intent(in) :: xx
	real(kind=r4),intent(in) :: x
	integer :: locate

	! Given an array xx(1:N), and given a value x, returns a value j such that x is between
	! xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
	! j = N is returned to indicate that x is out of range.

	integer :: n,jl,jm,Ju
	logical :: ascnd

	n=size(xx)
	ascnd = (xx(n) >= xx(1)) ! True if ascending order of table, false otherwise.
	jl=0                     ! Initialize Lower limit.
	ju=n+1                   ! Initialize upper limit.
	do
	  if(ju-jl <= 1) exit ! Repeat until this condition is satisfied
	  jm=(ju+jl)/2 ! compute a midpoint
	  if(ascnd .eqv. (x >= xx(jm))) then
	    jl=jm  
	  else
	    ju=jm
	  endif
	enddo
	if(x == xx(1)) then
	  locate=1
	else if (x == xx(n)) then
	  locate=n-1
	else
	  locate=jl
	endif
      end function locate


    end subroutine loglike_TL


