subroutine getRefractiveIndex2(lambda, grainType, mReal, mImg)
 use constants_mod
 use utils_mod
 use unix_mod

 implicit none
 real :: lambda
 real :: mReal, mImg
 integer, parameter :: maxpts = 100
 integer :: npts
 character(len=*) :: grainType
 real :: xArray(maxpts)
 real :: eps1para(maxpts)
 real :: eps2para(maxpts)
 real :: eps1perp(maxpts)
 real :: eps2perp(maxpts)
 real :: e1, e2, t
 real :: e1para, e1perp, e2para, e2perp
 integer :: i, j
 character(len=80) :: dataDirectory, filename

 call unixGetenv("TORUS_DATA",dataDirectory,i)
 dataDirectory = trim(dataDirectory)//"/"

 select case(grainType)

 case("graphite","draine")
    npts = 53
    filename = trim(dataDirectory)//"graphite_r.dat"
 open(20,file=filename, status="old", form="formatted")
 do i = 1, npts
  read(20,*) xArray(i), eps1para(i), eps2para(i), eps1perp(i), eps2perp(i)
 enddo
 close(20)

 xArray = xArray * 1.e-4
 
 if ((lambda*1.e-8 > xArray(1)).and.(lambda*1.e-8 < xArray(npts))) then
    call locate(xArray, npts, lambda*1.e-8, j)
    t = (lambda*1.e-8 - xArray(j))/(xArray(j+1)-xArray(j))
 else
    t = 0
    if (lambda*1.e-8 < xArray(1)) j = 1
    if (lambda*1.e-8 > xArray(npts)) j = npts-1
 endif

 e1para = eps1para(j) + t*(eps1para(j+1)-eps1para(j))
 e2para = eps2para(j) + t*(eps2para(j+1)-eps2para(j))
 

 e1perp = eps1perp(j) + t*(eps1perp(j+1)-eps1perp(j))
 e2perp = eps2perp(j) + t*(eps2perp(j+1)-eps2perp(j))
 
 e1 = 0.6666666*e1perp + 0.3333333*e1para
 e2 = 0.6666666*e2perp + 0.3333333*e2para

 mReal = oneByRootTwo * sqrt(sqrt(e1*e1 + e2*e2)+e1)
 mImg  = oneByRootTwo * sqrt(sqrt(e1*e1 + e2*e2)-e1)


case("silicate")
   npts = 30
    filename = trim(dataDirectory)//"silicate.dat"

 open(20,file=filename, status="old", form="formatted")
 do i = 1, npts
  read(20,*) xArray(i), eps1para(i), eps2para(i)
 enddo
 close(20)

 xArray = xArray * 1.e-4
 call locate(xArray, npts, lambda*1.e-8, j)
 t = (lambda*1.e-8 - xArray(j))/(xArray(j+1)-xArray(j))

 e1para = eps1para(j) + t*(eps1para(j+1)-eps1para(j))
 e2para = eps2para(j) + t*(eps2para(j+1)-eps2para(j))
 
 e1 = e1para
 e2 = e2para

 mReal = oneByRootTwo * sqrt(sqrt(e1*e1 + e2*e2)+e1)
 mImg  = oneByRootTwo * sqrt(sqrt(e1*e1 + e2*e2)-e1)

case DEFAULT

 write(*,'(a,a)') "! Grain type not recognised: ",grainType
 
end select

end subroutine getRefractiveIndex2
