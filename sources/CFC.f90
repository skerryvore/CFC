program NFCA
use iso_c_binding
use precision_kind
use CFC
use fwd
use utilities

implicit none
  include "mpif.h"
  integer,parameter :: LU_param=32
  integer,parameter :: LU_data=33
  integer,parameter :: LU_output=34
  integer,parameter :: LU_dummy=37
  integer,parameter :: nrow=9
  integer (kind = 4) :: nt2
  integer (kind = 4), dimension(:,:), allocatable :: ltri
  integer (kind = 4) :: lct(1)
  character(len = 150 ) :: data_file
  character(len = 200 ) :: param_file  
  character(len = 150 ) :: output_file
  character(len = 150 ) :: dem_file
  character(len = 150 ) :: filtered_file
  character(len = 200 ) :: filename
  character(len = 1024) :: line
  character(len = 20  ) :: dummy
  integer :: ndata
  integer :: icol
  integer :: i,j,io,kk,L,K,M, kk2
  real(kind=r8),allocatable,dimension(:,:) :: distances
  integer,allocatable,dimension(:,:) :: dem_topo,flt_topo
  real(kind=r8) :: rlon1,rlat1,rlon2, rlat2
  real(kind=r4)  :: range(2,1024)
  integer :: nd
  integer :: ierr, r1
  logical :: file_exists
  integer :: nmaps, ncsamp, ncloc
  real(kind=r4), dimension(:),allocatable :: maps 
  real(kind=r4) :: interp1
  real(kind=r4) :: dt
  real(kind=r4) :: mfitmin
  real(kind=r4) :: model_opt(50)
  real(kind=r4 ) :: cutwave
  integer :: mopt
  integer :: nproc
  integer :: iproc
  integer :: CIRCLE_DEF
  integer, parameter :: NCLOSEST=2, RADIUS =1, VORONOI=3 
  logical :: lroot  
  logical :: debug
  logical :: FIRST_PATH, BOREHOLE
  logical :: time_sequence(10)
  logical :: prntx
  real(kind=r4) :: t1,t2
  real(kind=r4) :: MTL, AFTA, FTLD(200), LKH 
  real(kind=r4) :: llTL, llFT 
  real(kind=r4), dimension(:), allocatable :: sorted_distances
  integer, dimension(:), allocatable :: idx
  integer :: sampleid
  CHARACTER(LEN=100), DIMENSION(200), TARGET :: stringArray
  TYPE(C_PTR), DIMENSION(200) :: stringPtrs

  ! List of variables used by TRIPACK
  integer(kind = 4), dimension(:), allocatable :: LIST, LPTR, LEND, NEAR, NEXT
  integer :: IER, LNEW
  real(kind=r8),dimension(:), allocatable :: DIST
  real(kind=8) :: wx1, wx2, wy1, wy2
  character(len = 80) :: title
  logical :: numbr
  integer (kind = 4), parameter :: lplt = 72
  real ( kind = 8 ), parameter :: pltsiz = 7.5D+00

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

  real(kind=c_double) :: dfGeoX, dfGeoY
  real(kind=c_double) :: elevation
  real(kind=c_double) :: rx, ry

  include "Cinterfaces.h"
  
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
  read(LU_param,*) filtered_file ! read name of dem file  
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

  ! open sample file list
  open(LU_data,file='./data/'//trim(adjustL(data_file)), status='unknown')
  
  ! read number of lines in the data file
  ndata=0
  read(LU_data,*) ! header
  do
    ndata = ndata + 1
    read(LU_data,*,iostat=io)
    if (io < 0) then
     ndata = ndata - 1
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
    sample(I)%s1name,&  
    sample(I)%lon,&
    sample(I)%lat,&  
    sample(I)%strat_old,&
    sample(I)%strat_young,&
    sample(I)%measured_elevation,&
    sample(I)%filepath
    sample(I)%id = I
  end do

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
  
  ! Get UTM Coordinates
  do I=1,ndata
    call LatLon2UTM(real(sample(I)%lon, kind=c_double),&
                    real(sample(I)%lat, kind=c_double),&
                    rx,ry)
    sample(I)%x = rx
    sample(I)%y = ry
  end do 

  do I=1,ndata
    call getElevation(real(sample(I)%lon,kind=c_double),&
                      real(sample(I)%lat,kind=c_double),&
                      "./DEM/"//trim(adjustL(dem_file))//C_NULL_CHAR,&
                      elevation) 
    sample(I)%dem_elevation = real(elevation)    
    call getElevation(real(sample(I)%lon,kind=c_double),&
                      real(sample(I)%lat,kind=c_double),&
                      "./DEM/"//trim(adjustL(filtered_file))//C_NULL_CHAR,&
                      elevation) 
    sample(I)%flt_elevation = real(elevation)
  enddo

  open(16, file="./DEM/Elevation_compare.txt",status='unknown')
  write(16,'(3a10,2a10,2a15)') "Measured ", "Extracted ", "Filtered ", "Lat ", "Lon ", "XUTM", "YUTM"
  do I = 1, ndata
    write(16,'(3i10,2f10.2, 2f15.2)') int(sample(I)%measured_elevation), sample(I)%dem_elevation, sample(I)%flt_elevation, &
              & sample(I)%lat, sample(I)%lon, sample(I)%x, sample(I)%y
  end do
  close(16)

  ! calculate distances. Distances are calculated in a horizontal plane.
  do I=1,ndata
    do J=1,ndata
      rlon1 = sample(I)%lon ; rlat1 = sample(I)%lat
      rlon2 = sample(J)%lon ; rlat2 = sample(J)%lat
      call GetDistanceBetweenPoints(rlon1, rlat1, rlon2, rlat2, distances(I,J))
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
          sample(I)%neighbours_offsets(kk) = sample(J)%dem_elevation-sample(J)%flt_elevation&
                                           - (sample(I)%dem_elevation - sample(I)%flt_elevation)
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
        sample(I)%neighbours_offsets(kk) = sample(idx(K))%dem_elevation-sample(idx(K))%flt_elevation&
                                         - (sample(I)%dem_elevation - sample(I)%flt_elevation)

        if(.not.BOREHOLE) ncsamp = ncsamp + 1
        K = K + 1

        deallocate(sorted_distances)
        deallocate(idx)
      end do

    end select
      
    end do
  end if

  if(CIRCLE_DEF == VORONOI) then 
    NCC = 0 
    ! Initialize arrays, LIST, LPTR, LEND, LNEW, NEAR, NEXT, DIST, IER
    IER = 0
    allocate(LIST(6 * ndata - 12))
    allocate(LPTR(6 * ndata - 12))
    allocate(LEND(ndata))
    allocate(NEAR(ndata))
    allocate(NEXT(ndata))
    allocate(DIST(ndata))

    ! Calculate delaunay triangulation using TRIPACK

    write(*, '(a)') "TRIPACK library"
    write(*, '(a,i6)') "The number of nodes is ", ndata

    ! Create the Delaunay triangulation (TRMESH), and test for errors
    call TRMESH(ndata, sample%x, sample%y,&
                LIST, LPTR, LEND,LNEW, NEAR, NEXT, DIST,ier) 
    if ( ier == -2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Error in TRMESH:'
      write ( *, '(a)' ) '  The first three nodes are collinear.'
      stop
    else if ( ier == -4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Error in TRMESH:'
      write ( *, '(a)' ) '  Invalid triangulation.'
      stop
    else if ( 0 < ier ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Error in TRMESH:'
      write ( *, '(a)' ) '  Duplicate nodes encountered.'
      stop
    else if ( ier > 0) then
      write ( *, '(a)' ) ' Error '
      stop 
    end if


    prntx = .true.
    allocate(ltri(nrow,2*ndata-5))
    call trprnt ( ncc, lcc, ndata, sample%x, sample%y, list, lptr, lend, prntx )
    call trlist ( ncc, lcc, ndata, list, lptr, lend, nrow, nt2, ltri, lct, ier )
    call trlprt ( ncc, lct, ndata, sample%x, sample%y, nrow, nt2, ltri, prntx )
    deallocate(ltri)

    ! Create an EPS file image of the triangulation
    sample%xnorm = sample%y / maxval(sample%y)
    sample%ynorm = sample%x / maxval(sample%x)
    wx1 = minval ( sample(1:ndata)%xnorm )
    wx2 = maxval ( sample(1:ndata)%xnorm )
    wy1 = minval ( sample(1:ndata)%ynorm )
    wy2 = maxval ( sample(1:ndata)%ynorm )

    !  Create an encapsulated PostScript file image of the triangulation.
    !
    open ( unit = lplt, file = 'tripack_prb.eps', status = 'unknown' )
    !
    !  NUMBR = TRUE iff nodal indexes are to be displayed.
    !
    numbr = ( ndata <= 200 )
    !
    !  Store a plot title.  It must be enclosed in parentheses.
    !
    title = '(Triangulation)'
   
    call trplot ( lplt, pltsiz, wx1, wx2, wy1, wy2, ncc, lcc, ndata, &
        sample%xnorm, sample%ynorm, LIST, LPTR, LEND, title, numbr, IER )
 
    close ( unit = lplt )
 
    if ( ier == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  TRPLOT has created a triangulation plot.'
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  TRPLOT failed with IER = ', ier
    end if

    ! Find neigbours of each sample
    allocate(NNABS(ndata))
    allocate(NPTR(ndata))
    allocate(NPTR1(ndata))
    allocate(NABOR(ndata))
    allocate(NBNOS(ndata))

    do I=1, ndata

      call find_node_neighbors(I,LIST, LPTR, LEND, NABOR, kk)
      sample(I)%nneighbours            = kk
      sample(I)%neighbours(1:kk)       = NABOR(1:kk)

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
        sample(I)%neighbours_offsets(L) = sample(J)%dem_elevation-sample(J)%flt_elevation&
                                        - (sample(I)%dem_elevation - sample(I)%flt_elevation)          
        print*, sample(I)%neighbours_offsets(L)
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
   
    ! Initialize all arrays (Must be done before doing calculation) 
   
    fwd_ages(1:kk)    = sample(I)%neighbours_ages(1:kk)
    fwd_ndata         = sample(I)%nneighbours
    fwd_MTL(1:kk)     = sample(I)%neighbours_MTL(1:kk)
    fwd_offsets(1:kk) = sample(I)%neighbours_offsets(1:kk)
    fwd_ncounts(1:kk) = sample(I)%neighbours_ncounts(1:kk)
    fwd_zeta(1:kk)    = sample(I)%neighbours_zeta(1:kk)
    fwd_rhodos(1:kk)  = sample(I)%neighbours_rhodos(1:kk)
    fwd_ntl(1:kk)     = sample(I)%neighbours_ntl(1:kk)
 
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
    range(2,4:(2*ttpoints-2)) = 150._r4 ! Temperature
    range(2,(2*ttpoints-1))   = 20._r4  ! Temperature
    
    mfitmin = 0.

    ! ATTENTION: NA calls the forward model which update some variables in the
    ! modules.
        
    call na(nd,range)

    write(*,*) "misfit:", mfitmin

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
          print*, AFTA, sample(I)%FTage
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
  write(LU_output,'(a15,11a10)') "Sample Name", "Lowest Misfit", "t1","t2","t3","t4","t5","T1","T2","T3","T4","T5" 
  do I=1,ndata  
    write(LU_output,'(a15,11f10.1)') sample(I)%s2name,sample(I)%Lwmisfit,sample(I)%bestpath(1,:),&
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
      if(K.gt.0) then
        sample(I)%Erate_rec(J) = (sample(I)%optimum_path(2,K+1) - sample(I)%optimum_path(2,K)) / &
    	  & (sample(I)%optimum_path(1,K+1) - sample(I)%optimum_path(1,K)) 
      end if
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
      &  real(sample(:)%x,c_double), real(sample(:)%y, c_double), int(ndata,i8))


    ! Create the .prj file associated to the shapefile
    ! The coordinate system is identical to the input topography file so
    ! we just need to copy the information from the topo.prj
    !open(71,file="./DEM/"//trim(adjustL(dem_file(1:index(dem_file,".asc",.false.)))//"prj"),status='old')
    !read(71,'(a)') line
    !close(71)
    !open(71,file="Temp"//trim(adjustL(filename))//".prj",status="Unknown")
    !write(71,'(a)') trim(adjustL(line))
    !close(71)

    ! Create a dbf file
    do I = 1, ndata
      stringArray(I) = sample(I)%s1name//C_NULL_CHAR
      stringPtrs(I) = C_LOC(stringArray(I))
    end do

    call createDBF("Temp"//trim(adjustL(filename))//".dbf"//achar(0),&
      & sample(:)%id, int(ndata,i8),  real(sample(:)%x, c_double),&
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
  integer :: I,J,K
  integer,parameter :: NPOINTS=4
  integer :: nd
  real(kind=r4) :: TIME(4),TEMP(4),ORDERED_TIME_SEQ(4)
  real(kind=r4) :: ftage(500),ftldistrib(500,200),ftldmean(500)  
  real(kind=r4) :: NA_param(nd)
  real(kind=r4) :: misfit
  real(kind=r4) :: LKH_FTAGE(NSAMPLEMAX)
  real(kind=r4) :: LKH_TLDistrib(NSAMPLEMAX) 
  real(kind=r4) :: LKH_sample(NSAMPLEMAX)
  real(kind=r4) :: total_LKH

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
      misfit=1e22
      return
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
  misfit=abs(total_LKH)        

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
      use CFC

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
      integer :: idx

      real(kind=r4) :: models (nd,*)
      real(kind=r4) :: misfit (ntot)
      character*(*) header
      character*5 :: tempo


      write(tempo,'(i5)') model_id
      open (87,file='./results/'//trim(adjustL(run_name))//'/NA_results_sample'//trim(adjustL(tempo))//'.txt',status='unknown')
      write (87,*) "Total number of samples:", sample(model_id)%nneighbours + 1
      write (87,'(80a)') ("-", I = 1, 80)
      write (87,300) "Name", "lat", "lon","Elev.meas", "Elev.dem","Elev.flt","FTage", "MTL"
      write (87,'(80a)') ("-", I = 1, 80)
      write (87,301)   sample(model_id)%s1name,&
                     sample(model_id)%lat,&
                     sample(model_id)%lon,&
                     int(sample(model_id)%measured_elevation),&
                     sample(model_id)%dem_elevation,&
                     sample(model_id)%flt_elevation,&
                     sample(model_id)%FTage,&
                     sample(model_id)%MTL

      do i=1,sample(model_id)%nneighbours
     
        idx = sample(model_id)%neighbours(i)
        write (87,301) sample(idx)%s1name,&
                     sample(idx)%lat,&
                     sample(idx)%lon,&
                     int(sample(idx)%measured_elevation),&
                     sample(idx)%dem_elevation,&
                     sample(idx)%flt_elevation,&
                     sample(idx)%FTage,&
                     sample(idx)%MTL
      end do


      write (87,'(80a)') ("-", I = 1, 80)
      write (87, *) "Model with lowest misfit"
      write (87, 302) sample(model_id)%optimum_LKH, sample(model_id)%optimum_path(1,1:3), sample(model_id)%optimum_path(2,:)   
      write (87,'(80a)') ("-", I = 1, 80)
      write (87, *) "List of models" 
      write (87,'(80a)') ("-", I = 1, 80)
     
      write (87, 300) "Misfit", "t1", "t2", "t3", "T1", "T2", "T3", "T4" 
      write (87,'(80a)') ("-", I = 1, 80)
      do i=1,ntot
      
        write (87,302) misfit(i),(models(k,i),k=1,nd)
      
      enddo
      close (87)
      
      300 format(8a10)
      301 format(a10, 2f10.1, 3i10, 2f10.1)
      302 format(E15.8,100f10.1)

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


