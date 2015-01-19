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
  
  interface 
    subroutine getElevation(dfGeox, dfGeoY, Filename, elevation) bind(C, name="getElevation")
      use iso_c_binding
      implicit none
      real(kind=c_double ) :: dfGeoX, dfGeoY
      character(kind=c_char) :: Filename(*)
      real(kind=c_double) :: elevation
    end subroutine
  end interface
 
  interface 
    subroutine GetDistanceBetweenPoints(X1, Y1, X2, Y2, distance) bind(C, name="GetDistanceBetweenPoints")
      use iso_c_binding
      implicit none
      real(kind=c_double ) :: X1,Y1,X2,Y2
      real(kind=c_double) :: distance
    end subroutine
  end interface
  
  interface 
    subroutine LatLon2UTM(X1, Y1, X2, Y2) bind(C, name="LatLon2UTM")
      use iso_c_binding
      implicit none
      real(kind=c_double ) :: X1,Y1,X2,Y2
    end subroutine
  end interface
