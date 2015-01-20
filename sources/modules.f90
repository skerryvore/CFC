module precision_kind
  implicit none
  integer, parameter   :: i4=SELECTED_INT_KIND(4)
  integer, parameter   :: i8=SELECTED_INT_KIND(8)
  integer, parameter   :: r4=SELECTED_REAL_KIND(6,37)
  integer, parameter   :: r8=SELECTED_REAL_KIND(15,307)

end module

module CFC
  use precision_kind
  
  real(kind=r4) :: circle
  integer,parameter :: ttpoints=4 
  type database
    sequence 
    integer :: id
    character :: s1name*50
    character :: s2name*50
    character :: filepath*200
    real(kind=r8) :: x
    real(kind=r8) :: y
    real(kind=r8) :: xnorm
    real(kind=r8) :: ynorm
    integer :: z
    real(kind=r8) :: lat
    real(kind=r8) :: lon
    character :: country*50
    integer :: NCOUNTS
    integer :: NTL
    real(kind=r4) :: zeta
    real(kind=r4) :: rhod
    real(kind=r4) :: nd
    real(kind=r4) :: FTage
    real(kind=r4) :: FTage_err
    real(kind=r4) :: FTLD(17) 
    real(kind=r4) :: MTL
    real(kind=r4) :: MTL_Err
    real(kind=r4) :: MTL_std
    real(kind=r4) :: stdev
    real(kind=r4) :: stdev_err
    real(kind=r4) :: offset
    integer :: NS(100)
    integer :: NI(100)
    real(kind=r4) :: TL(400)
    integer :: neighbours(500) 
    integer :: nneighbours
    real(kind=r4) :: neighbours_ages(500)
    real(kind=r4) :: neighbours_ages_err(500)
    real(kind=r4) :: neighbours_offsets(500)
    integer       :: neighbours_ncounts(500)
    integer       :: neighbours_ntl(500)
    real(kind=r4) :: neighbours_zeta(500)
    real(kind=r4) :: neighbours_rhodos(500)
    integer       :: neighbours_NS(500,100)
    integer       :: neighbours_NI(500,100)
    real(kind=r4) :: neighbours_MTL(500)
    real(kind=r4) :: neighbours_MTL_err(500)
    real(kind=r4) :: neighbours_TL(500,400)
    real(kind=r4) :: neighbours_FTLD(500,17)
    real(kind=r4) :: bestgeotherm
    real(kind=r4) :: bestmod_geotherm
    integer :: family(500)
    integer :: nfamily
    real(kind=r4) :: misfit
    real(kind=r4) :: Lwmisfit
    real(kind=r4) :: best_mod(2,4) 
    real(kind=r4) :: bestpath(2,4)
    real(kind=r4) :: temp_at_time(500,2)
    real(kind=r4) :: strat_young
    real(kind=r4) :: strat_old
    real(kind=r4) :: measured_elevation  
    integer :: dem_elevation
    integer :: flt_elevation
    real(kind=r4) :: time_rec(20) 
    real(kind=r4) :: temp_rec(20)
    real(kind=r4) :: Erate_rec(20)
    real(kind=r4) :: optimum_path(2,4)
    real(kind=r4) :: optimum_LKH
    real(kind=r4) :: optimum_geotherm
  end type database
  type (database),dimension(:),allocatable :: sample
end module CFC

module fwd
use precision_kind
integer,parameter :: NSAMPLEMAX=100
integer,parameter :: NCOUNTMAX=100
integer,parameter :: NTLMAX=400
real(kind=r4)  :: fwd_offsets(NSAMPLEMAX+1)
integer  :: fwd_ndata
integer  :: fwd_NS(NSAMPLEMAX,NCOUNTMAX)
integer  :: fwd_NI(NSAMPLEMAX,NCOUNTMAX)
real(kind=r4) :: fwd_TL(NSAMPLEMAX,NTLMAX)
real(kind=r4) :: fwd_ages(NSAMPLEMAX)
real(kind=r4) :: fwd_ages_error(NSAMPLEMAX)
real(kind=r4) :: fwd_MTL(NSAMPLEMAX)
real(kind=r4) :: fwd_MTL_error(NSAMPLEMAX)
integer :: fwd_ncounts(NSAMPLEMAX)
real(kind=r4) :: fwd_zeta(NSAMPLEMAX)
real(kind=r4) :: fwd_rhodos(NSAMPLEMAX)
integer :: fwd_ntl(NSAMPLEMAX)

real(kind=r4) :: path(2,4)
integer :: model_id
logical :: search_geo
real(kind=r4) :: geotherm
character*150 :: run_name

end module

! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

module qsort_c_module
use precision_kind

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A)
  real(kind=r4), intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real(kind=r4), intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real(kind=r4) :: temp
  real(kind=r4) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module

module SPindex 
  
! From Numerical recipes in Fortran 90
! The overloading has been removed and only the required routines are
! present.

INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: SP = KIND(1.0)

contains

SUBROUTINE indexx(arr,index)
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: arr
  INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
  INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
  ! Indexes an array arr, i.e., outputs the array index of length N such that arr(index(j ))
  ! is in ascending order for j = 1, 2, . . . , N . The input quantity arr is not changed.
  REAL(SP) :: a
  INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
  INTEGER(I4B), DIMENSION(NSTACK) :: istack
  n=assert_eq(size(index),size(arr),'indexx_sp')
  index=arth(1,1,n)
  jstack=0
  l=1
  r=n
  do
    if (r-l < NN) then
      do j=l+1,r
        indext=index(j)
        a=arr(indext)
        do i=j-1,l,-1
          if (arr(index(i)) <= a) exit
          index(i+1)=index(i)
        end do
        index(i+1)=indext
      end do
      if (jstack == 0) RETURN
      r=istack(jstack)
      l=istack(jstack-1)
      jstack=jstack-2
    else
      k=(l+r)/2
      call swap(index(k),index(l+1))
      call icomp_xchg(index(l),index(r))
      call icomp_xchg(index(l+1),index(r))
      call icomp_xchg(index(l),index(l+1))
      i=l+1
      j=r
      indext=index(l+1)
      a=arr(indext)
      do
        do
          i=i+1
          if (arr(index(i)) >= a) exit
        end do
        do
          j=j-1
          if (arr(index(j)) <= a) exit
        end do
        if (j < i) exit
        call swap(index(i),index(j))
      end do
      index(l+1)=index(j)
      index(j)=indext
      jstack=jstack+2
      if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
      if (r-i+1 >= j-l) then
        istack(jstack)=r
        istack(jstack-1)=i
        r=j-1
      else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
      end if
    end if
  end do

CONTAINS

SUBROUTINE icomp_xchg(i,j)
  INTEGER(I4B), INTENT(INOUT) :: i,j
  INTEGER(I4B) :: swp
  if (arr(j) < arr(i)) then
    swp=i
    i=j
    j=swp
  end if
END SUBROUTINE icomp_xchg

END SUBROUTINE indexx

FUNCTION arth(first,increment,n)
  INTEGER(I4B), INTENT(IN) :: first,increment,n
  INTEGER(I4B), DIMENSION(n) :: arth
  INTEGER(I4B) :: k,k2,temp
  INTEGER(I4B),PARAMETER :: NPAR_ARTH=16, NPAR2_ARTH=8
  if (n > 0) arth(1)=first
  if (n <= NPAR_ARTH) then
    do k=2,n
      arth(k)=arth(k-1)+increment
    end do
  else
    do k=2,NPAR2_ARTH
      arth(k)=arth(k-1)+increment
    end do
    temp=increment*NPAR2_ARTH
    k=NPAR2_ARTH
    do
      if (k >= n) exit
      k2=k+k
      arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
      temp=temp+temp
      k=k2
    end do
  end if
END FUNCTION arth

FUNCTION assert_eq(n1,n2,string)
  !  Report and die if integers not all equal (used for size checking).
  CHARACTER(LEN=*), INTENT(IN) :: string
  INTEGER, INTENT(IN) :: n1,n2
  INTEGER :: assert_eq
  if (n1 == n2) then
    assert_eq=n1
  else
    write (*,*) 'nrerror: an assert_eq failed with this tag:', &
      string
    STOP 'program terminated by assert_eq'
  end if
END FUNCTION assert_eq

SUBROUTINE nrerror(string)
  !  Report a message, then die.
  CHARACTER(LEN=*), INTENT(IN) :: string
  write (*,*) 'nrerror: ',string
  STOP 'program terminated by nrerror'
END SUBROUTINE nrerror

SUBROUTINE swap(a,b)
  ! Swap the contents of a and b.
  INTEGER(I4B), INTENT(INOUT) :: a,b
  INTEGER(I4B) :: dum
  dum=a
  a=b
  b=dum
END SUBROUTINE swap

end module SPindex


module utilities
  use precision_kind

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

end module utilities



