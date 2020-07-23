Program Main

!gfortran -o xdis_PW distance_PW.f90 tecio64.a -lstdc++

  implicit none

  integer, parameter :: im = 351, jm = 121, km = 121, nzone = 2

  integer, dimension(nzone) :: img, jmg, kmg

  integer ::i, j, k, nn, nz, kbb, ii, imin, imax

  integer ::  iob1, iob2, job1, job2, cont, ntotal ! limites del obstaculo y cantidad total de nodos

  real, dimension(im, jm, km, nzone) :: x, y, z, wd, wd_sa ! wd_sa: distancia al borde solido mas cercano
  real::  wda1, wda2, wda3, one, cdesd, H
  real::  uta, ren, duma, dumb 
  real:: d1, d2, d3, d4, d5, d6 ! Variables para comparar distancias a los distintos bordes solidos

! variables to enable writing of TecPlot binary (*.plt) files

  integer (kind = 4)           :: TecIni, TecDat, TecZne, TecNod, TecFil
  integer (kind = 4)           :: TecEnd
  integer (kind = 4)           :: VIsDouble = 0
  integer (kind = 4)           :: Debug = 1
  
  integer (kind = 4)           :: III
  character (len = 256) :: filename
  character (len = 1)          :: nullchr = char(0)



  one = 1.0E+00
  H = 0.377643505

  print*, 'Lectura de ind3dmg...'
  !---------------------
  open (unit = 21, file = 'ind3dmg.dat', form = 'formatted')
  read (unit=21,fmt=*) duma, dumb
  do nz = 1, nzone
     read (unit = 21, fmt = *) img(nz), jmg(nz), kmg(nz)
     print*, nz, img(nz), jmg(nz), kmg(nz)
  end do
  close(unit = 21)
 

  !print*, 'Lectura del archivo nodos_obstaculo...'
  !open(unit = 41, file = "nodos_obstaculo", form = "formatted")
  !  read(41, fmt = "(2I5)") iob1, iob2
  !  read(41, fmt = "(2I5)") job1, job2
  !close(unit = 41)
  !print*, 'nodos obstaculo: '
  !print*, 'iob1: ', iob1
  !print*, 'iob2: ', iob2
  !print*, 'job1: ', job1
  !print*, 'job2: ', job2
  !print*, ' '  
!

!****************************************************************
! Lectura de la malla
!****************************************************************	

open(1,file="grid.dat",form="unformatted")
  
  do nn=1,nzone
     read(1) nz
   print *,nz
     read(1) img(nz), jmg(nz), kmg(nz)
   print*, img(nz), jmg(nz), kmg(nz)
     read(1) (((x(i,j,k,nz), i = 1,img(nz)), j = 1,jmg(nz)), k = 1,kmg(nz))
     read(1) (((y(i,j,k,nz), i = 1,img(nz)), j = 1,jmg(nz)), k = 1,kmg(nz))
     read(1) (((z(i,j,k,nz), i = 1,img(nz)), j = 1,jmg(nz)), k = 1,kmg(nz))
  end do
  close(1)

!****************************************************************
! Calculo de la distancia a la pared segun la zona
!**************************************************************** 

! Convencion: x direccion del flujo e y transversal

! nzob1: SW del obstaculo
! nzob2: W del obstaculo
! nzob3: NW del obstaculo
! nzob4: S del obstaculo
! nzob5: OBSTACULO
! nzob6: N del obstaculo
! nzob7: SO del obstaculo
! nzob8: O del obstaculo
! nzob9: NO del obstaculo

wd = 0.0

!----------------------------------------------
nz=1
cdesd = 0.65

!******************************************************************
! CALCULO DE LA DISTANCIA A BORDE SOLIDO
!*****************************************************************

cont = 0

i=1
do nz=1,nzone
do k=1,kmg(nz)
  do j=1,jmg(nz)

      d1=99999.9 ! jm wall - main channel
      d2=99999.9 ! bottom - main channel
      d3=99999.9 ! interface wall - main channel
      d4=99999.9 ! interface corner
      d5=99999.9 ! bottom - floodplain
      d6=99999.9 ! wall - floodplain

      if(nz.eq.1) then
      
        d1 = sqrt((x(i,j,k,nz)-x(i,jmg(1),k,1))**2+(y(i,j,k,nz)-y(i,jmg(1),k,1))**2+(z(i,j,k,nz)-z(i,jmg(1),k,1))**2)
        d2 = sqrt((x(i,j,k,nz)-x(i,j,1,1))**2+(y(i,j,k,nz)-y(i,j,1,1))**2+(z(i,j,k,nz)-z(i,j,1,1))**2)

        if (z(i,j,k,nz).lt.H) then
          d3 = sqrt((x(i,j,k,nz)-x(i,1,k,1))**2+(y(i,j,k,nz)-y(i,1,k,1))**2+(z(i,j,k,nz)-z(i,1,k,1))**2)
        end if
      
      end if


      d4 = sqrt((x(i,j,k,nz)-x(i,1,51,1))**2+(y(i,j,k,nz)-y(i,1,51,1))**2+(z(i,j,k,nz)-z(i,1,51,1))**2)

      if(nz.eq.2) then
        d5 = sqrt((x(i,j,k,nz)-x(i,j,1,nz))**2+(y(i,j,k,nz)-y(i,j,1,nz))**2+(z(i,j,k,nz)-z(i,j,1,nz))**2)
        d6 = sqrt((x(i,j,k,nz)-x(i,1,k,nz))**2+(y(i,j,k,nz)-y(i,1,k,nz))**2+(z(i,j,k,nz)-z(i,1,k,nz))**2)
      end if

      wd_sa(:,j,k,nz) = min(d1,d2,d3,d4,d5,d6)

  end do
end do
end do


ntotal = im*jm*km

print*, 'Malla: ',ntotal
print*, 'Nodos recorridos: ', cont
print*, ' '

!******************************************************************
! CALCULO DEL VOLUMEN DE ELEMENTO PARA LES
!*****************************************************************

  do nz = 1, nzone
     do k = 1,kmg(nz)
        do j = 1,jmg(nz)
           do i = 1,img(nz)

              if(i<img(nz)) then
                wda1 = sqrt((x(i+1,j,k,nz)-x(i,j,k,nz))**2+(y(i+1,j,k,nz)-y(i,j,k,nz))**2+(z(i+1,j,k,nz)-z(i,j,k,nz))**2)
              else
                wda1 = sqrt((x(i,j,k,nz)-x(i-1,j,k,nz))**2+(y(i,j,k,nz)-y(i-1,j,k,nz))**2+(z(i,j,k,nz)-z(i-1,j,k,nz))**2)
              end if

              if(j<jmg(nz)) then
                wda2 = sqrt((x(i,j+1,k,nz)-x(i,j,k,nz))**2+(y(i,j+1,k,nz)-y(i,j,k,nz))**2+(z(i,j+1,k,nz)-z(i,j,k,nz))**2)
              else
                wda2 = sqrt((x(i,j-1,k,nz)-x(i,j,k,nz))**2+(y(i,j-1,k,nz)-y(i,j,k,nz))**2+(z(i,j-1,k,nz)-z(i,j,k,nz))**2)
              end if

              if(k<kmg(nz))then
                wda3 = sqrt((x(i,j,k+1,nz)-x(i,j,k,nz))**2+(y(i,j,k+1,nz)-y(i,j,k,nz))**2+(z(i,j,k+1,nz)-z(i,j,k,nz))**2)
              else
                wda3 = sqrt((x(i,j,k,nz)-x(i,j,k-1,nz))**2+(y(i,j,k,nz)-y(i,j,k-1,nz))**2+(z(i,j,k,nz)-z(i,j,k-1,nz))**2)
              end if

              cdesd = 0.65*max(wda1,wda2,wda3)
              wd(i,j,k,nz) = min(wd_sa(i,j,k,nz),cdesd)

              !!!!!!!!!!!!!!!!!!!!!!
              ! FOR URANS:
              !!!!!!!!!!!!!!!!!!!!!!
              wd(i,j,k,nz) = wd_sa(i,j,k,nz) 

          end do
        end do
      end do
    end do

   
  ! Escribien archivo grid

   print *, 'Comenzando a escribir archivo grid'
   print *, '...'
   open  (unit = 41, file = 'grid', form='unformatted')
   do nz = 1, nzone

      write  (unit = 41) (((x (i,j,k,nz),i=1,img(nz)), &
                                       j=1,jmg(nz)),k=1,kmg(nz))
      write  (unit = 41) (((y (i,j,k,nz),i=1,img(nz)), &
                                       j=1,jmg(nz)),k=1,kmg(nz))
      write  (unit = 41) (((z (i,j,k,nz),i=1,img(nz)), &
                                       j=1,jmg(nz)),k=1,kmg(nz))
      write  (unit = 41) (((wd(i,j,k,nz),i=1,img(nz)), &
                                       j=1,jmg(nz)),k=1,kmg(nz))
   end do
   close(unit = 41)

  print *,'Archivo grid escrito!'
  print *, '...'



print *, 'Comenzando a escribir archivo TECPLOT'
print *, '...'

!------------TECPLOT BEGIN
write(filename, fmt = '(a,i6.6)')'distance'
	
    I = TecIni('Turbine'//NULLCHR,	&   ! title of file
	'X, Y, Z, WD'//NULLCHR,	&   ! list of variables
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
    ! write each variable
    ! the last argument in the following calls indicates
    ! the precision of the the variable,
    ! 0 = single, 1 = double
		I   = TecDat(III,x(1:img(nz),1:jmg(nz),1:kmg(nz),nz),0)
		I   = TecDat(III,y(1:img(nz),1:jmg(nz),1:kmg(nz),nz),0)
		I   = TecDat(III,z(1:img(nz),1:jmg(nz),1:kmg(nz),nz),0)
    I   = TecDat(III,wd(1:img(nz),1:jmg(nz),1:kmg(nz),nz),0)
end do
I   = TecEnd()



print *, 'Archivo distance.plt escrito!'
print *, '...'
!------------TECPLOT END

  print *,'end!'

!--------------------------------------------

END Program Main
