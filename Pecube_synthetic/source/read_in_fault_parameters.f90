!----------------------------

subroutine read_in_fault_parameters (fault,nfault,xlon1,xlat1,xlon2,xlat2,xl,yl,zl,timeend, &
                                     nd,range,param)

! this subroutines reads the fault information from the fault_parameters.txt file
! and generates the necessary complementary information to be used for fault
! kinematics and movement

use definitions

implicit none

integer nfault,i,j,k,icol,jcol,kcol,jd
type(faulttype) fault(nfault)
double precision xn,yn,xyn,x1,y1,x2,y2,timeend
double precision xlon1,xlat1,xlon2,xlat2,xl,yl,zl
character line*1024

real*4 range(2,*),param(*)
integer nd

if (nfault.eq.0) return

      open (76,file='input/fault_parameters.txt',status='old')
      open (77,status='scratch')
    1 read (76,'(a1024)',end=2) line
        if (line(1:1).ne.'$'.and. line(1:1).ne.' ') then
          if (scan(line,'$').ne.0) then
            do i=scan(line,'$'),1024
            line(i:i)=' '
            enddo
          endif
        k=1
          do j=1,1024
            if (line(j:j).eq.' '.or.line(j:j).eq.',') then
            if (j.ne.k) write (77,'(a)') line(k:j-1)
            k=j+1
            endif
          enddo
        endif
      goto 1
    2 close (76)
      rewind (77)

read (77,*)
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
nd=nd+1
read (line(1:icol-1),*) range(1,nd)
read (line(icol+1:1024),*) range(2,nd)
x1=param(nd)
elseif (jcol.ne.0) then
read (line(jcol+1:1024),*) jd
x1=param(jd)
else
read (line,*) x1
endif
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
nd=nd+1 
read (line(1:icol-1),*) range(1,nd)
read (line(icol+1:1024),*) range(2,nd)
y1=param(nd)
elseif (jcol.ne.0) then
read (line(jcol+1:1024),*) jd
y1=param(jd)
else    
read (line,*) y1
endif
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
nd=nd+1 
read (line(1:icol-1),*) range(1,nd)
read (line(icol+1:1024),*) range(2,nd)
x2=param(nd)
elseif (jcol.ne.0) then
read (line(jcol+1:1024),*) jd
x2=param(jd)
else    
read (line,*) x2
endif
read (77,'(a1024)') line
icol=scan(line,':')
jcol=scan(line,'#')
if (icol.ne.0) then
nd=nd+1 
read (line(1:icol-1),*) range(1,nd)
read (line(icol+1:1024),*) range(2,nd)
y2=param(nd)
elseif (jcol.ne.0) then
read (line(jcol+1:1024),*) jd
y2=param(jd)
else    
read (line,*) y2
endif
x1=(x1-xlon1)/(xlon2-xlon1)*xl
y1=(y1-xlat1)/(xlat2-xlat1)*yl
x2=(x2-xlon1)/(xlon2-xlon1)*xl
y2=(y2-xlat1)/(xlat2-xlat1)*yl
  do i=1,nfault
  read (77,*) fault(i)%n
  if (fault(i)%n.lt.0) then
  allocate (fault(i)%x(4),fault(i)%y(4))
    do k=1,4
    read (77,'(a1024)') line
    icol=scan(line,':')
    jcol=scan(line,'#')
    if (icol.ne.0) then
    nd=nd+1
    read (line(1:icol-1),*) range(1,nd)
    read (line(icol+1:1024),*) range(2,nd)
    fault(i)%x(k)=param(nd)
    elseif (jcol.ne.0) then
    read (line(jcol+1:1024),*) jd
    fault(i)%x(k)=param(jd)
    else
    read (line,*) fault(i)%x(k)
    endif
    enddo
  else
  allocate (fault(i)%x(fault(i)%n),fault(i)%y(fault(i)%n))
    do k=1,fault(i)%n
    read (77,'(a1024)') line
    icol=scan(line,':')
    jcol=scan(line,'#')
    if (icol.ne.0) then
    nd=nd+1 
    read (line(1:icol-1),*) range(1,nd)
    read (line(icol+1:1024),*) range(2,nd)
    fault(i)%x(k)=param(nd)
    elseif (jcol.ne.0) then
    read (line(jcol+1:1024),*) jd
    fault(i)%x(k)=param(jd)
    else    
    read (line,*) fault(i)%x(k)
    endif
    read (77,'(a1024)') line
    icol=scan(line,':')
    jcol=scan(line,'#')
    if (icol.ne.0) then
    nd=nd+1 
    read (line(1:icol-1),*) range(1,nd)
    read (line(icol+1:1024),*) range(2,nd)
    fault(i)%y(k)=param(nd)
    elseif (jcol.ne.0) then
    read (line(jcol+1:1024),*) jd
    fault(i)%y(k)=param(jd)
    else    
    read (line,*) fault(i)%y(k)
    endif
    fault(i)%y(k)=fault(i)%y(k)+zl
    enddo
  fault(i)%x1=x1;fault(i)%y1=y1;fault(i)%x2=x2;fault(i)%y2=y2
  endif
  read (77,*) fault(i)%nstep
  allocate (fault(i)%timestart(fault(i)%nstep),fault(i)%timeend(fault(i)%nstep), &
            fault(i)%velo(fault(i)%nstep))
    do k=1,fault(i)%nstep
    read (77,'(a1024)') line
    icol=scan(line,':')
    jcol=scan(line,'#')
    kcol=scan(line,'*')
    if (icol.ne.0) then
    nd=nd+1 
    read (line(1:icol-1),*) range(1,nd)
    read (line(icol+1:1024),*) range(2,nd)
    fault(i)%timestart(k)=param(nd)
    elseif (jcol.ne.0) then
    read (line(jcol+1:1024),*) jd
    fault(i)%timestart(k)=param(jd)
    elseif (kcol.ne.0) then
    if (k.eq.1) stop 'cannot use wild card for first time step'
    !fault(i)%timestart(k)=fault(i)%timeend(k-1)
    fault(i)%timestart(k)=timeend-fault(i)%timeend(k-1) ! Mod Thibaud
    else
    read (line,*) fault(i)%timestart(k)
    endif
    read (77,'(a1024)') line
    icol=scan(line,':')
    jcol=scan(line,'#')
    if (icol.ne.0) then
    nd=nd+1 
    read (line(1:icol-1),*) range(1,nd)
    read (line(icol+1:1024),*) range(2,nd)
    fault(i)%timeend(k)=param(nd)
    elseif (jcol.ne.0) then
    read (line(jcol+1:1024),*) jd
    fault(i)%timeend(k)=param(jd)
    else    
    read (line,*) fault(i)%timeend(k)
    endif
    read (77,'(a1024)') line
    icol=scan(line,':')
    jcol=scan(line,'#')
    if (icol.ne.0) then
    nd=nd+1 
    read (line(1:icol-1),*) range(1,nd)
    read (line(icol+1:1024),*) range(2,nd)
    fault(i)%velo(k)=param(nd)
    elseif (jcol.ne.0) then
    read (line(jcol+1:1024),*) jd
    fault(i)%velo(k)=param(jd)
    else    
    read (line,*) fault(i)%velo(k)
    endif
    fault(i)%timestart(k)=timeend-fault(i)%timestart(k)
    fault(i)%timeend(k)=timeend-fault(i)%timeend(k)
    enddo
  enddo
close (77)
 
  do i=1,nfault
  if (fault(i)%n.gt.0) &
  allocate (fault(i)%xs(fault(i)%n-1),fault(i)%ys(fault(i)%n-1))
  enddo

call calculate_fault_parameters (fault,nfault)
  
  do i=1,nfault
  if (fault(i)%n.gt.0) then
  if (fault(i)%x(fault(i)%n).gt.fault(i)%x(1) .and. &
      fault(i)%y(fault(i)%n).lt.fault(i)%y(1)) &
          fault(i)%velo(1:fault(i)%nstep)=-fault(i)%velo(1:fault(i)%nstep)
  if (fault(i)%x(fault(i)%n).lt.fault(i)%x(1) .and. &
      fault(i)%y(fault(i)%n).gt.fault(i)%y(1)) &
          fault(i)%velo(1:fault(i)%nstep)=-fault(i)%velo(1:fault(i)%nstep)
  endif
  enddo
  
  do i=1,nfault
  xn=fault(i)%y2-fault(i)%y1
  yn=fault(i)%x1-fault(i)%x2
  xyn=sqrt(xn**2+yn**2)
  fault(i)%xn=xn/xyn
  fault(i)%yn=yn/xyn
  enddo

return

end subroutine read_in_fault_parameters
