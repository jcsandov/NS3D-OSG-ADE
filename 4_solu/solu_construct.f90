Program Main

  ! gfortran -o mksolu solu_construct.f90
  
  implicit none

  integer, parameter :: im = 351, jm = 121, km = 121, nzone =2

  integer, dimension(im, jm, km) :: ln

  integer, dimension(nzone) :: img, jmg, kmg

  integer ::i, j, k, nn, nz, kbb, m , l, tmp_size, me

  integer ::  iob1, iob2, job1, job2, cont, ntotal ! limites del obstaculo y cantidad total de nodos

  real, dimension(im, jm, km, nzone) :: x, y, z, wd
  real, dimension(im, jm, km, nzone) :: p,u,v,w,q5,xnut

  real, dimension(im, jm, km, nzone) :: xdum, ydum, zdum

  real, dimension(:,:), allocatable :: qtmp, qntmp
  real, dimension(:), allocatable :: nutmp
  
  integer :: dum1, dum2

me = 5
p=0.0E+00
u=1.0E+00
v=0.0E+00
w=0.0E+00
q5=0.0E+00
xnut=0.005E+00

  print*, 'Lectura del archivo ind3dmg...'
  print*, ' '
  open(1, file = "../3_distance/ind3dmg.dat", form = "formatted")
  read(1,*) dum1, dum2
  do nz = 1,nzone
    read(1,*) img(nz), jmg(nz),kmg(nz)
    print *, img(nz), jmg(nz), kmg(nz)
  end do
  close(1)

  print*, 'Lectura del archivo grid...'
  print*, ' '
  open(1,file="../3_distance/grid",form="unformatted")
   do nz = 1, nzone
     read(1) (((x(i,j,k,nz), i = 1,img(nz)), j = 1,jmg(nz)), k = 1,kmg(nz))
     read(1) (((y(i,j,k,nz), i = 1,img(nz)), j = 1,jmg(nz)), k = 1,kmg(nz))
     read(1) (((z(i,j,k,nz), i = 1,img(nz)), j = 1,jmg(nz)), k = 1,kmg(nz))
     read(1) (((wd(i,j,k,nz), i = 1,img(nz)), j = 1,jmg(nz)), k = 1,kmg(nz))
     end do
  close(1)

!  print*, 'Lectura del archivo nodos_obstaculo...'
!  open(unit = 41, file = "../nodos_obstaculo", form = "formatted")
!    read(41, fmt = "(2I5)") iob1, iob2
!    read(41, fmt = "(2I5)") job1, job2
!  close(unit = 41)
!  print*, 'nodos obstaculo: '
!  print*, 'iob1: ', iob1
!  print*, 'iob2: ', iob2
!  print*, 'job1: ', job1
!  print*, 'job2: ', job2
!  print*, ' '  



do nz = 1, nzone
do k= 1, kmg(nz)
  do j = 1,jmg(nz)
    do i = 1,img(nz)
      ! SOLID BOUNDARY CONDITION

      if(nz==2) u(i,j,k,nz) = 0.8
      if(wd(i,j,k,nz)==0) u(i,j,k,nz) = 0.0
 !     if(i.ge.iob1.and.i.le.iob2.and.j.ge.job1.and.j.le.job2) u(i,j,k,nz) = 0.0
      ! zero velocity at the air layer
!      if(z(i,j,k,nz).gt.1.0) u(i,j,k,nz) = 0.0

    end do
  end do
end do    
end do

!-------------------------------------------------------


          open  (unit = 42, file = 'solu', form = 'unformatted')
          do nz = 1, nzone

             ! temporary arrays 
             ! 
             tmp_size = img(nz) * jmg(nz) * kmg(nz)
       
             allocate(qtmp(1:me,tmp_size), &
                     qntmp(1:me,tmp_size) )
             ! allocate temp array for turbulence viscosity
             ! 
             allocate(nutmp(tmp_size))

       qtmp = 0.0E+00
       qntmp = 0.0E+00
       nutmp = 0.0E+00

       do k = 1, kmg(nz); do j = 1, jmg(nz); do i = 1, img(nz)
        ln(i,j,k) = (k - 1) * img(nz) * jmg(nz) &
              + (j - 1) * img(nz)      &
              +  i 
       end do; end do; end do

        do i = 1,img(nz)
        do j = 1,jmg(nz)
        do k = 1,kmg(nz)
          l = ln(i,j,k)
          qtmp(1,l)=p(i,j,k,nz)
          qtmp(2,l)=u(i,j,k,nz)
          qtmp(3,l)=v(i,j,k,nz)
          qtmp(4,l)=w(i,j,k,nz)
          qtmp(5,l)=q5(i,j,k,nz)
          nutmp(l)=xnut(i,j,k,nz)
        end do
        end do
        end do

       qntmp(1:me,1:tmp_size) = qtmp(1:me,1:tmp_size)

             do m = 1, me
                write  (unit = 42)  (qtmp(m,l), l = 1, tmp_size)
             end do
                write  (unit = 42)  (nutmp(l),  l = 1, tmp_size)
             do m = 1, me
                write  (unit = 42)  (qntmp(m,l),l = 1, tmp_size)
             end do

             deallocate (qtmp, qntmp, nutmp)
          
      end do
          close (unit = 42)



  print *,'end!'



END Program Main
