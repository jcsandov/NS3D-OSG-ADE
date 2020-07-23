Program Main

!gfortran -o xgrid get_grid.f90 tecio64.a -lstdc++

  implicit none

  integer, parameter :: im = 351, jm = 121, km = 121, nzone = 2

  !integer, dimension(im, jm, km, nzone) :: solid

  integer, dimension(nzone) :: img, jmg, kmg

  integer ::i, j, k, nn, nz, kbb, ii, imin, imax

  double precision, dimension(im, jm, km, nzone) :: xx, yy, zz, xxdum
  real, dimension(im, jm, km, nzone) :: x, y, z, xdum
  real::  one, cdesd, wd_sa
  real::  uta, ren, duma, dumb 

! variables to enable writing of TecPlot binary (*.plt) files

  integer (kind = 4)           :: TecIni, TecDat, TecZne, TecNod, TecFil
  integer (kind = 4)           :: TecEnd
  integer (kind = 4)           :: VIsDouble = 0
  integer (kind = 4)           :: Debug = 1
  
  integer (kind = 4)           :: III
  character (len = 256) :: filename
  character (len = 1)          :: nullchr = char(0)



  one = 1.0E+00

  !---------------------
  open (unit = 21, file = 'ind3dmg.dat', form = 'formatted')
  read (unit=21,fmt=*) duma, dumb
  do nz = 1, nzone
     read (unit = 21, fmt = *) img(nz), jmg(nz), kmg(nz)
     print*, nz, img(nz), jmg(nz), kmg(nz)
  end do
  close(unit = 21)




	
open(2,file="../0_grid_PW/CaseIV2.dat")
do nn=1,6
read(2,*)
end do

do nz=1,nzone
print*,'nz',nz
	do nn=1,5
	read(2,*)
	end do
!do j=1,jmg(nz)  ! Reads every cross-stream plane
!	do i=img(nz),1,-1; do k=1,kmg(nz)	
do k=1,kmg(nz)
print*,k
  do j=1,jmg(nz); do i=1,img(nz)
		read(2,*) xx(i,j,k,nz),yy(i,j,k,nz),zz(i,j,k,nz)
	end do; end do

end do
end do !nz

close(2)


x=real(xx)
y=real(yy)
z=real(zz)


print*, 'preparandose para escribir grid.dat'

open(1,file="grid.dat",form="unformatted") 
  do nn=1,nzone
     write(1) nn
print *,nn
     write(1) img(nn), jmg(nn), kmg(nn)
 print*, img(nn), jmg(nn), kmg(nn)
     write(1) (((x(i,j,k,nn), i = 1,img(nn)), j = 1,jmg(nn)), k = 1,kmg(nn))
     write(1) (((y(i,j,k,nn), i = 1,img(nn)), j = 1,jmg(nn)), k = 1,kmg(nn))
     write(1) (((z(i,j,k,nn), i = 1,img(nn)), j = 1,jmg(nn)), k = 1,kmg(nn))
  end do
  close(1)


!------------TECPLOT BEGIN
write(filename, fmt = '(a,i6.6)')'grid_organized'
	
    I = TecIni('Erosion'//NULLCHR,	&   ! title of file
	'X, Y, Z'//NULLCHR,	&   ! list of variables
	trim(filename)//'.plt'//NULLCHR,&   ! output file name
	'.'//NULLCHR,					&
	Debug,							&
	VIsDouble)
do nz=1,nzone

		I = TecZne('Zone'//NULLCHR,    &     
			img(nz),			           &
			jmg(nz),			           &
			kmg(nz),			           &
			'BLOCK'//NULLCHR,		   &
			NULLCHR//NULLCHR)
			! total number of points
		III = img(nz) * jmg(nz) * kmg(nz)
		I   = TecDat(III,xx(1:img(nz),1:jmg(nz),1:kmg(nz),nz),1)
		I   = TecDat(III,yy(1:img(nz),1:jmg(nz),1:kmg(nz),nz),1)
		I   = TecDat(III,zz(1:img(nz),1:jmg(nz),1:kmg(nz),nz),1)
end do

    I   = TecEnd()

!------------TECPLOT END


  
 !     open  (unit = 41, file = 'grid', form='unformatted')
 !     do nz = 1, nzone
 !        write  (unit = 41) (((x (i,j,k,nz),i=1,img(nz)), &
 !                                         j=1,jmg(nz)),k=1,kmg(nz))
 !        write  (unit = 41) (((y (i,j,k,nz),i=1,img(nz)), &
 !                                         j=1,jmg(nz)),k=1,kmg(nz))
 !        write  (unit = 41) (((z (i,j,k,nz),i=1,img(nz)), &
 !                                         j=1,jmg(nz)),k=1,kmg(nz))
 !         write  (unit = 41) (((wd(i,j,k,nz),i=1,img(nz)), &
 !                                          j=1,jmg(nz)),k=1,kmg(nz))
 ! 	end do
  	close(unit = 41)


  print *,'end!'

!--------------------------------------------

END Program Main
