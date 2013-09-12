subroutine read_mtx(input_file,ncounts,ntl,zeta,rhod,nd,FTage,FTage_error,&
                    &MTL,MTL_error,STDEV,STDEV_error,NS,NI,TL,sample_name)
use precision_kind

implicit none

integer   :: LU_mtx
integer   :: nconstraints
integer   :: ncounts
integer   :: ntl
real(kind=r4)   :: nd
integer   :: I
real(kind=r4) :: rhod
real(kind=r4) :: zeta
real(kind=r4) :: FTage,FTage_error
real(kind=r4) :: MTL,MTL_error
real(kind=r4) :: STDEV,STDEV_error
integer :: NS(100),NI(100)
real(kind=r4) :: TL(200)
character :: sample_name*50
character :: input_file*200

LU_mtx=89

open(LU_mtx,file=trim(adjustL(input_file)),status="unknown")
read(LU_mtx,*) sample_name
read(LU_mtx,*)
read(LU_mtx,*) nconstraints,ntl,ncounts,zeta,rhod,nd
do I=1,nconstraints
read(LU_mtx,*)
enddo

read(LU_mtx,*) FTage,FTage_error
read(LU_mtx,*) MTL,MTL_error
read(LU_mtx,*) STDEV,STDEV_error

do I=1,ncounts
  read(LU_mtx,*) NS(I), NI(I)
enddo


do I=1,ntl
  read(LU_mtx,*) TL(I)
enddo

end subroutine read_mtx
