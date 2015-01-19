
module precision_kind

  integer, parameter   :: i4=SELECTED_INT_KIND(4)
  integer, parameter   :: i8=SELECTED_INT_KIND(8)
  integer, parameter   :: sp=SELECTED_REAL_KIND(6,37)
  integer, parameter   :: dp=SELECTED_REAL_KIND(15,307)
end module precision_kind

subroutine generate_data(rho, age, zeta, pdf, vals, sampleid, mtl)
  use precision_kind
  implicit none

  real(kind = dp) :: rho, age

  integer, dimension(:), allocatable :: Ns, Ni 
  integer :: NsNi, Nc, totco, I, nloop, nconstraints
  integer :: sampleid
  real(kind = sp) :: zeta
  real(kind = dp) :: ld = 1.555125e-4
  real(kind = sp) :: rhod
  real(kind = sp) :: prob
  real(kind = dp) :: centage, relerr, std_err
  integer(kind = i8) :: SEED, SEED2
  real(kind = sp) :: dum
  real(kind = sp), dimension(17) :: pdf, vals
  real(kind = sp), dimension(:), allocatable :: TL 
  character(len=100) :: str
  real(kind = sp) mtl, mtl_std


  ! Random number seed
  call system_clock(SEED)
  call random_seed(SEED)
  
  ! Number of crystals
  Nc = 30

  allocate(Ni(Nc))
  allocate(Ns(Nc))

  ! Specify rho_d and Nd for dosimeter... from age equation
  ! The factor 2.0e6 is for the geometry factor and the decay constant
  ! in the correct units (as used to calculate the central age)
  rhod = 2.0e6 / zeta / rho / ld * (exp(ld * age) - 1.0)
  totco = 2000
  ! probability in binomial distribution
  prob = rho / (1. + rho)

  ! For Nc crystals, generate synthetic Ns and Ni count data using binomial 
  ! distribution, conditional on total counts Ns + Ni
  do I = 1, Nc
    NsNi = 0
    do while (NsNi == 0)
      NsNi = int(ran3(SEED)*1000.0)
    end do
    Ns(I) = 0
    do while (Ns(I) == 0)
      Ns(I) = int(bnldev(prob, NsNi, SEED))
      if(Ns(I) == 0) print*, Ns(I), NsNi
      Ni(I) = NsNi - Ns(I)
    end do
  end do

  ! Calculate central age
  centage = 0.0
  call centralage(centage, relerr, std_err, Nc, Ns, Ni, zeta, totco, rhod)

  ! Pick 100 random lengths
  allocate(TL(100))

  call random_from_distrib(pdf, vals, TL, 100)
  
  ! Write mtx files
  nconstraints = 0 ! Not used, should be 0

  write(str,'(i3)') sampleid
  open(66, file="RUN00/SyntheticData/Sample"//trim(adjustL(str))//".mtx", status="unknown")
  write(*,*) "Sample"//trim(adjustL(str))//" Central age = ", centage, " Input age = ", age

  ! Write the name of the sample
  write(66,*) "Sample"//trim(adjustL(str))
  write(66,'(i4)') -999
  write(66,'(3i4,f7.1,f10.1,i5)') nconstraints, size(TL),&
              size(Ns), zeta, rhod,&
              totco 

  do I = 1, nconstraints
    write(66,*) -999 ! not used
  end do

  ! FT age, Ft age error
  write(66,*) age, 0.05*age
  mtl = sum(TL) / 100
  mtl_std = sqrt(sum(TL**2)/100 - mtl**2)
  write(66,*) mtl, 0.05*mtl
  write(66,*) mtl_std, mtl_std * 0.05 

  do I = 1, Nc
    write(66, *) Ns(I), Ni(I)
  end do

  do I = 1, 100
   write(66,*) TL(I)
  end do
  
  close(66)

!  ! Write QTQT file
!  open(66, file="RUN00/SyntheticData/Sample"//trim(adjustL(str))//".qtqt", status="unknown")
!
! !Line 1
! write(66,'(a)') "Sample"//trim(adjustL(str))
! 
! !Line 2
! LAT = 0.0 ; LON = 0.0; ELEVATION = 0.0
! write(string(1),'(f10.2)') LAT
! write(string(2),'(f10.2)') LON
! write(string(3),'(f10.0)') ELEVATION
! write(66,'(a,1x,a,1x,a)') (trim(adjustL(string(I))),I=1,3)
! 
!
! close(66)
  deallocate(Ni)
  deallocate(Ns)
  deallocate(TL)

contains 

subroutine centralAge(centage, relerr, std_err, Nc, Ns, Ni, zeta, totco, rho_d)
  use precision_kind
  implicit none

  integer :: Nc
  integer :: Ns(Nc), Ni(Nc)
  integer :: totco
  real(kind = dp) :: centage, relerr, std_err
  real(kind = sp) :: zeta	
  real(kind = sp) :: rho_d
  integer :: I, J, nstot, nitot
  real(kind = dp) :: zbar, zsd, sigma, eta, mu, sumwy, sumw, ssqwy, ssq, eta2, semu

  real(kind = dp), dimension(:), allocatable :: Z, Y, w
  integer, dimension(:), allocatable :: m

  if( Nc == 1) then
    i = 1
    centage = 1. / (1.55125e-4)*log(1. + 1.55125e-10 * 0.5 * zeta * rho_d * ( Ns(i) + 0.5 ) / ( Ni(i) + 0.5 ))
    relerr  = sqrt(1. / ( Ns(i) + 0.5 ) + 1. / ( Ni(i) + 0.5 ) + 1. / ( totco + 0.5))
    std_err = centage * relerr
    return
  end if

  zbar = 0.0
  zsd = 0.0
  nstot = 0
  nitot = 0
  sigma = 0.0

  allocate(Z(Nc))
  allocate(Y(Nc))
  allocate(w(Nc))
  allocate(m(Nc))

  do I = 1, Nc
    m(I) = Ns(I) + Ni(I)
    y(I) = Ns(I) / real(m(I), dp)
    Z(I) = log((Ns(I) + 0.5 ) / ( Ni(I) + 0.5))
    zbar = zbar + Z(I)
    zsd = zsd + Z(I)**2 
    nstot = nstot + Ns(I)
    nitot = nitot + Ni(I)   
  end do

  zsd = sqrt(zsd - zbar**2 / Nc) / (Nc - 1)
  zbar = zbar / Nc
  eta = nstot / real( nstot + nitot, dp)
  mu = log(real(nstot/nitot, dp))
  sigma = 0.6 * zsd

  do I = 1, 20
    sumw = 0.0
    sumwy = 0.0
    ssq = (sigma * eta * ( 1. - eta))**2
    do J = 1, Nc
      w(J) = m(J) / ( eta * (1.0 - eta ) + ( m(J) - 1 ) * ssq)
      sumwy = sumwy + w(J) * y(J)
      sumw = sumw + w(J)
    end do

    eta = sumwy / sumw
    mu = log(eta/(1.0 - eta))
    ssqwy = 0.0
    do J = 1, Nc
      ssqwy = ssqwy + (w(J)*(y(J) - eta))**2
      sigma = sigma * sqrt(ssqwy / sumw)	
    end do
  end do

  eta2 = (eta * (1.0 - eta))**2
  semu = sqrt(1. / (sumw * eta2) + 1. / totco)
  centage = 1. / (1.55125e-4) * log(1.+1.55125e-10*0.5*zeta*rho_d*exp(mu))
  relerr = 100. * sigma
  std_err = centage * semu

  deallocate(Z)
  deallocate(Y)
  deallocate(w)
  deallocate(m)
end subroutine centralAge

function bnldev(pp, n, idum)
  use precision_kind
  implicit none
  real(kind = sp) :: pp
  integer :: n
  integer(kind = i8) :: idum 
  real(kind = sp), parameter :: PI=3.141592654
  integer :: J
  integer, save :: nold = -1
  real(kind = sp ) ::am, em, g, angle, p, bnldev, sq, t, y
  real(kind = sp ), save :: pold = -1
  real(kind = sp ), save :: pc, plog, pclog, en, oldg

  p = merge(pp, 1.0_sp - pp, pp <= 0.5_sp)

  am = n*p
  if (n < 25) then
    bnldev = 0.
    do J = 1, n
      if( ran3(idum) < p ) bnldev = bnldev + 1.0
    end do
  else if (am < 1.0) then
    g = exp(-am)
    t = 1.0
    do J = 0, n
      t = t * ran3(idum)
      if(t < g) exit
    end do
    bnldev = merge(j, n, j <= n)
  else
    if ( n .ne. nold) then
      en = n
      oldg = gammln(real(en + 1.0_sp, dp))
      nold = n
    end if
    if ( p .ne. pold) then
      pc = 1.0_sp - p
      plog = log(p)
      pclog = log(pc)
      pold = p
    end if
    sq = sqrt(2.0_sp*am*pc)
    do 
      angle = PI * ran3(idum)
      y = tan(angle)
      em = sq*y+am
      if(em < 0.0 .or. em >= en+1.0_sp) cycle
      em = int(em)
      t = 1.2_sp * sq * (1.0_sp + y**2) * exp(oldg - gammln(real(em + 1.0_sp, dp)) - &
        &  gammln(real(en - em + 1.0_sp, dp)) + em * plog + (en -em) * pclog)
      if (ran3(idum) <= t ) exit
    end do
    bnldev = em
  end if

  if (p .ne. pp ) bnldev = n - bnldev
end function bnldev

function gammln(xx)
  use precision_kind
  implicit none
  ! returns the log natural of the gamma function. From numerical recipes.
  real(kind = dp) :: gammln
  real(kind = dp) :: xx
  integer :: J
  real(kind = dp) :: ser, tmp, x, y
  real(kind = dp), parameter  :: stp = 2.5066282746310005D0 
  real(kind = dp), dimension(6), parameter  :: cof = (/ 76.180091732947146d0, -86.50532032941677d0, &
    & 24.01409824083091d0, -1.231739572450155d0, &
    & 0.12086509738661793d-2,-0.5395239384953d-5 /) 
  x = xx
  y = x
  tmp = x + 5.5d0
  tmp = (x + 0.5d0) * log(tmp) - tmp
  ser = 1.000000000190015d0
  do J = 1, 6
    y = y + 1.d0
    ser = ser + cof(J) / y
  end do  
  gammln = tmp + log(stp * ser / x)
end function gammln


function ran3(idum)
  implicit none
  integer :: idum
  real ran3
  integer, parameter :: MBIG = 1000000000, MSEED = 161803398, MZ = 0 
  real, parameter :: FAC = 1./MBIG
  integer :: I, iff, ii, inext, inextp, k
  integer :: mj, mk, ma(55)
  save iff, inext, inextp, ma
  data iff /0/

  if(idum.lt.0.or.iff.eq.0) then
    iff = 1
    mj=abs(MSEED - abs(idum))
    mj=mod(mj,MBIG)
    ma(55)=mj
    mk=1
    do I=1, 54
      ii = mod(21*I, 55)
      ma(ii) = mk
      mk = mj - mk
      if(mk.lt.MZ) mk = mk + MBIG
      mj = ma(ii)
    end do
    do k = 1,4
      do i = 1,55
	ma(i) = ma(i) - ma(1+mod(i+30,55))
	if(ma(i).lt.MZ) ma(i) = ma(i) + MBIG
      enddo
    end do
    inext = 0
    inextp = 31
    idum = 1
  end if


  inext = inext + 1
  if(inext.eq.56) inext = 1
  inextp=inextp + 1
  if(inextp.eq.56) inextp = 1
  mj = ma(inext) - ma(inextp)
  if(mj.lt.MZ) mj = mj + MBIG
  ma(inext) = mj
  ran3 = mj*FAC
end function ran3

!-----------------------------------------------------------------------
!
! subroutine random_from_distrib
!
! This is a simple routine to extract a given number N value X from a given
! probability density function (pdf)
!
! X can take discrete values (x(i), i=1, i=lenght(pdf)
! and we know from the pdf the probability p(i) that x takes the value x(i),
! such that Sum(p(i), i=1, lenght(pdf)) = 1.0
! 
subroutine random_from_distrib(pdf, vals, X, N)
  use precision_kind
  implicit none

  integer :: I, J, N
  real(kind = sp) :: X(N), R
  real(kind = sp), dimension(:) :: pdf, vals 
  real(kind = sp) ,dimension(:), allocatable :: cdf 

  ! First check that the values in the pdf array sum to 1
  ! For some reason that I have not figured out yet the sum of the pdf
  ! distribution out of ketcham routine is not 1. It is one only for a number of
  ! bins equal to 20 but Pecube uses 17 bins. The values needs to be normalized
  pdf = pdf / sum(pdf)

  ! Check that pdf and vals have the same size
  if(size(pdf) /= size(vals)) STOP "PDF and VALS don't have the same size"

  ! Calculate the cumulative density function cdf.
  allocate(cdf(size(pdf)))

  do I = 1, size(pdf)
    CDF(i) = sum(pdf(1:i))
  end do

  ! Generate a random number from a uniform distribution between 0 and 1
  ! Here I use the random number generator associated to the compiler but
  ! one can do better if needed.

  do I = 1, N
    call random_number(R)
    X(I) = pdf(1)
    do J = 2, size(cdf)
      if(CDF(j-1).lt.R.and.CDF(j).ge.R) X(I) = vals(J)
    end do
  end do

end subroutine random_from_distrib


end subroutine
