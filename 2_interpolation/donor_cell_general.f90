PROGRAM MAIN
! The donor-cell search code
! Input: grid.dat & bcs.dat
! 	bcs.dat: The control file which supplies the boundary conditions
!	The orgnization of the files is as following:
!-----------------------------------------------------------------------------
!	Number_of_blocks
!	For block=1,number_of_blocks
!		Number_of_interfaced_blocks,which_block
!		boundary conditions for six surfaces
!		For interface, 0 means interface, 
!		1 means outer boundary
!	END
!		e.g.
!			2, 2,3
!			0 0 1 1 1 1
!		The example here means for this block has two interfaced blocks
!		block #2 & block #3. 
!		The i=1 and i=imax surface are interface which the interface
!		  information need be searched.
!
!	The second part of the bcs.dat is the blanking information
!	For block=1,number_of_blocks
!		number_of_contained_blocks,contained_block_#
!		containing_surface_condition in six directions
!			(if ==0, then search, otherwise, skip)
!		is,ie,js,je,ks,ke of the blanked area
!	END
!		Here, if number_of_contained_blocks==0, then no blanking

!----------------------------------------------------------------------------
!	grid.dat is orgnized as
!	For block=1,number_of_blocks
!		imax,jmax,kmax
!		(((x(i,j,k),i=1,imax),j=1,jmax),k=1,kmax)
!	END

  !USE linkedlist3d
  !USE kdtree2_module

  implicit none
  integer::i1,i2,i3,nzone,nmax,MAXITS,nzones
  parameter(i1=351,i2=101,i3=121,nzones=2,nmax=2000000,MAXITS=200)
  real::dis(6),rst(3), t1, t2

  real,dimension(:,:,:,:),allocatable::x,y,z
  integer,dimension(:,:,:,:),allocatable::searchflag,iblank

  logical::check
  real::xmin(nzones),ymin(nzones),xmax(nzones),ymax(nzones),zmin(nzones),zmax(nzones)

  integer,dimension(:,:),allocatable::hostst,nodest,layer
  real,dimension(:,:),allocatable::coeft
  integer:: iii, n_int_total, n_int_max

  integer,dimension(nzones)::nint,zint(nzones,nzones),nblank,zblank(nzones,nzones)
  integer,dimension(6,nzones)::bctype,blktype
  integer,dimension(nzones)::isblk,ieblk,jsblk,jeblk,ksblk,keblk,n_index_nz

  integer,dimension(6)::istart,iend,jstart,jend,kstart,kend
  integer,dimension(nzones)::im,jm,km,nii(i1,i2,i3,nzones,nzones),njj(i1,i2,i3,nzones,nzones)
  integer::nz,m,i,j,k,iface,nn,ni,ii,jj,kk,inside,n_index,nl,its,itss,kkk
  real::minx,miny,minz,maxx,maxy,maxz
  real::point(3),cell(8,3)
  real::p(3),v000(3),v100(3),v010(3),v001(3),v101(3),v011(3),v110(3),v111(3)
  common/intp/v000,v100,v010,v001,v101,v011,	&
       &   v110,v111,p

! allocatable local variables for kdtree search:
!type (kdtree2), pointer    :: tree
!type (kdtree2_result), allocatable :: results(:)
!real  , dimension(:,:), allocatable :: xt
integer :: ibck, ifwd, jbck, jfwd, kbck, kfwd, lt
integer :: ibase, jbase, kbase

!allocate(results(2))


  nii=0
  njj=0
  print *,'prepare to read bcs.dat'

  call cpu_time(t1)

  open(1,file='bcs.dat')
  read(1,*) nzone
  if(nzone>nzones)then
     print *,'The nzone is greater than the maximum zone number in this code!'
     print *,'Please use a larger number of nzones!'
  end if
  allocate(x(i1,i2,i3,nzone),y(i1,i2,i3,nzone),z(i1,i2,i3,nzone),	&
       & searchflag(i1,i2,i3,nzone),iblank(i1,i2,i3,nzone),		&
       & hostst(4,nmax),coeft(3,nmax),layer(nmax,nzone),	&
       & nodest(4,nmax))

  layer=0
  searchflag=0


  do nz=1,nzone
     read(1,*) nint(nz),(zint(nz,m),m=1,nint(nz))
     read(1,*) (bctype(iface,nz),iface=1,6)
  end do
  do nz=1,nzone
     read(1,*) nblank(nz),(zblank(nz,m),m=1,nblank(nz))
     if(nblank(nz).ne.0)then
        read(1,*) (blktype(iface,nz),iface=1,6)
        read(1,*) isblk(nz),ieblk(nz),jsblk(nz),jeblk(nz),ksblk(nz),keblk(nz)
     end if
  end do

  close(1)
print *,'finish reading bcs.dat'
  open(1,file='grid.dat',form='unformatted')
  do nn=1,nzone
     read(1) nz
     read(1) im(nz),jm(nz),km(nz)
     if(im(nz)>i1)then
        print *,'The number of grid nodes along i direction is too large!'
        print *,'Please change i1 value in this code and recompile!'
        stop
     end if
     if(jm(nz)>i2)then
        print *,'The number of grid nodes along j direction is too large!'
        print *,'Please change i2 value in this code and recompile!'
        stop
     end if
     if(km(nz)>i3)then
        print *,'The number of grid nodes along k direction is too large!'
        print *,'Please change i3 value in this code and recompile!'
        stop
     end if
     read(1) (((x(i,j,k,nz),i=1,im(nz)),j=1,jm(nz)),k=1,km(nz))
     read(1) (((y(i,j,k,nz),i=1,im(nz)),j=1,jm(nz)),k=1,km(nz))
     read(1) (((z(i,j,k,nz),i=1,im(nz)),j=1,jm(nz)),k=1,km(nz))

     xmin(nz)=minval(x(1:im(nz),1:jm(nz),1:km(nz),nz))
     xmax(nz)=maxval(x(1:im(nz),1:jm(nz),1:km(nz),nz))
     ymin(nz)=minval(y(1:im(nz),1:jm(nz),1:km(nz),nz))
     ymax(nz)=maxval(y(1:im(nz),1:jm(nz),1:km(nz),nz))
     zmin(nz)=minval(z(1:im(nz),1:jm(nz),1:km(nz),nz))
     zmax(nz)=maxval(z(1:im(nz),1:jm(nz),1:km(nz),nz))

  end do
  close(1)

  open(1,file='interface.dat',form='unformatted')
  iblank=0
  do nz=nzone,1,-1
     if(nblank(nz).ne.0)then
        do k=ksblk(nz),keblk(nz)-1
        do j=jsblk(nz),jeblk(nz)-1
        do i=isblk(nz),ieblk(nz)-1
           iblank(i,j,k,nz)=1
        end do
        end do
        end do
     end if
  end do
  
iii = 0

  do nz=1,nzone
     print *,'nz',nz,im(nz),jm(nz),km(nz)

   kstart(1)=2
     kend(1)=km(nz)-1
     jstart(1)=1
     jend(1)=jm(nz)
     istart(1)=1
     iend(1)=3

     kstart(2)=2
     kend(2)=km(nz)-1
     jstart(2)=1
     jend(2)=jm(nz)
     istart(2)=im(nz)-2
     iend(2)=im(nz)


     kstart(3)=2
     kend(3)=km(nz)-1
     jstart(3)=1
     jend(3)=1
     istart(3)=1
     iend(3)=im(nz)

     kstart(4)=2
     kend(4)=km(nz)-1
     jstart(4)=jm(nz)
     jend(4)=jm(nz)
     istart(4)=1
     iend(4)=im(nz)

      if (nz.eq.1) then
              kstart(3)=52
              kend(3)=120
      end if

     kstart(5)=1
     kend(5)=2
     jstart(5)=1
     jend(5)=jm(nz)
     istart(5)=1
     iend(5)=im(nz)

     kstart(6)=km(nz)-1
     kend(6)=km(nz)
     jstart(6)=1
     jend(6)=jm(nz)
     istart(6)=1
     iend(6)=im(nz)



     ! The following block search and calculate the interpolation coefficients
     ! of the block interface. The program use iblanking array to decide 
     ! whether a cell is a blanked cell. Searchflag is the array stores the
     ! information whether the point is searched or not.
     ! The program call subroutine search to decide the relationship
     ! between a point and cell. If the returned variable inside==1, then
     ! the point is within the cell and subroutine newt is called to decide
     ! the interpolation coefficients.
     ! The searching point is first transfered to point, and the donor cell
     ! geometries are transfered to cell(8,3). When the point is inside
     ! the cell, the coefficients are stored. If a point lies in different
     ! blocks, all of them are searched and stored. The corresponding array
     ! to store the information is nlayer. nlayer==the number of donor blocks.
     n_index=0
     do iface=1,6
	     if(bctype(iface,nz).eq.0.and.nint(nz).ge.1) then

                 searchloopo:do nn=1,nint(nz)
                    ni=zint(nz,nn)


           !do i=istart(iface),iend(iface)
           do k=kstart(iface),kend(iface)
           do j=jstart(iface),jend(iface)
           do i=istart(iface),iend(iface)

              if(searchflag(i,j,k,nz).ne.1) then
                 !searchflag(i,j,k,nz)=1
                 point(1)=x(i,j,k,nz)
                 point(2)=y(i,j,k,nz)
                 point(3)=z(i,j,k,nz)
          
                ibck = 1 ;      ifwd = im(ni)-1
                jbck = 1;       jfwd = jm(ni)-1
                kbck = 1;       kfwd = km(ni)-1

                if (ni.eq.1) then
                        kbck = 51
                        kfwd = 120
                end if

                       searchloop:do kk=kbck,kfwd 
                       do jj=jbck,jfwd 
                       do ii=ibck,ifwd 
                          inside=0
                          if(iblank(ii,jj,kk,ni).eq.0)then
                             cell(1,1)=x(ii,jj,kk,ni)
                             cell(1,2)=y(ii,jj,kk,ni)
                             cell(1,3)=z(ii,jj,kk,ni)
                             cell(2,1)=x(ii+1,jj,kk,ni)
                             cell(2,2)=y(ii+1,jj,kk,ni)
                             cell(2,3)=z(ii+1,jj,kk,ni)
                             cell(3,1)=x(ii+1,jj+1,kk,ni)
                             cell(3,2)=y(ii+1,jj+1,kk,ni)
                             cell(3,3)=z(ii+1,jj+1,kk,ni)
                             cell(4,1)=x(ii,jj+1,kk,ni)
                             cell(4,2)=y(ii,jj+1,kk,ni)
                             cell(4,3)=z(ii,jj+1,kk,ni)
                             cell(5,1)=x(ii,jj,kk+1,ni)
                             cell(5,2)=y(ii,jj,kk+1,ni)
                             cell(5,3)=z(ii,jj,kk+1,ni)
                             cell(6,1)=x(ii+1,jj,kk+1,ni)
                             cell(6,2)=y(ii+1,jj,kk+1,ni)
                             cell(6,3)=z(ii+1,jj,kk+1,ni)
                             cell(7,1)=x(ii+1,jj+1,kk+1,ni)
                             cell(7,2)=y(ii+1,jj+1,kk+1,ni)
                             cell(7,3)=z(ii+1,jj+1,kk+1,ni)
                             cell(8,1)=x(ii,jj+1,kk+1,ni)
                             cell(8,2)=y(ii,jj+1,kk+1,ni)
                             cell(8,3)=z(ii,jj+1,kk+1,ni)
                             call search(point,cell,inside,dis)

                             if(inside.eq.1) then
                               print *,'================================================================'
                 	             print *,'succes for',i,j,k,nz
                               print *, '---------------------------------------------------------------'
                               print *,'donor',ii,jj,kk,ni
                               print *,' '

                                iii = iii + 1
                                searchflag(i,j,k,nz)=1
                                n_index=n_index+1

                                nodest(1,iii)=i
                                nodest(2,iii)=j
                                nodest(3,iii)=k
                                nodest(4,iii)=nz

                                layer(n_index,nz)=layer(n_index,nz)+1
                                nl=layer(n_index,nz)

                                hostst(1,iii)=ii
                                hostst(2,iii)=jj
                                hostst(3,iii)=kk
                                hostst(4,iii)=ni

                                do m=1,3
					                     !!print *,dis(2*m-1),dis(2*m)
                                   rst(m)=dis(2*m-1)/(dis(2*m-1)+dis(2*m))
                                end do
                                
                                !print *,'--------------------------------------------------------------'
                                !print *,'rst(m) inicial: ', rst(1), rst(2), rst(3)
                                !print *, '   '
                                
                                do m=1,3
                                   v000(m)=cell(1,m)
                                   v100(m)=cell(2,m)
                                   v010(m)=cell(4,m)
                                   v001(m)=cell(5,m)
                                   v101(m)=cell(6,m)
                                   v011(m)=cell(8,m)
                                   v110(m)=cell(3,m)
                                   v111(m)=cell(7,m)
                                   p(m)=point(m)
                                end do

                                !call newt(rst,3,check,MAXITS,its,itss)
                                !print *,'--------------------------------------------------------------'
                                !print *,'rst(m) despues de newt: ', rst(1), rst(2), rst(3)
                                !print *, '   '

                                if(its<MAXITS)then
                                   do m=1,3
                                      coeft(m,iii)=rst(m)
                                   end do
                                else
                                   coeft(m,iii)=dis(2*m-1)/(dis(2*m-1)+dis(2*m))
                                end if ! its
                                print *,'--------------------------------------------------------------'
                                print *,'Coeficientes calculados'
                                print *,coeft(1:3,iii)
					                      print *,' '



                                exit searchloop
                             end if ! inside== 1
                          end if ! iblank == 0
                       end do
                       end do
                       end do searchloop


              end if !searchflag
              if(searchflag(i,j,k,nz).ne.1.and.nn.eq.nint(nz))then
                 print *,'Search failed at point',i,j,k,nz
                 
                 stop
              end if


           end do ! i
           end do ! j 
           end do ! k


                 end do searchloopo  !do nn=1,nint(nz)

        end if    ! bctype
     end do       ! iface
     


 

n_index_nz(nz) = n_index

  end do !NZ


  n_int_total = 0
  n_int_max = 0
  do nz = 1,nzone
	n_int_total = n_int_total + n_index_nz(nz)
	n_int_max = max(n_int_max,n_index_nz(nz))
  end do

  print *,n_int_total,n_int_max
  if (n_int_total>nmax) print*,'ERROR n_int_total>nmax'

  write (unit = 1) n_int_total, n_int_max

       iii = 0
       do nz = 1, nzone
	      print *,'ZONE nz=',nz,'n_index_nz=',n_index_nz(nz)
          write (unit = 1) n_index_nz(nz)
       do i = 1, n_index_nz(nz)
          iii = iii + 1

          write (unit = 1) nodest(1,iii),nodest(2,iii),nodest(3,iii),nodest(4,iii)
          write (unit = 1) hostst(1,iii),hostst(2,iii),hostst(3,iii),hostst(4,iii)
          write (unit = 1)  coeft(1,iii), coeft(2,iii), coeft(3,iii)
       end do
       end do


  call cpu_time(t2)
  print *, 'time', t2 - t1

  close(1)

  open(2,file='interface.check',form='formatted')
  write (2,*) n_int_total, n_int_max

  iii = 0
  do nz = 1, nzone
     write (2,*) n_index_nz(nz)
     do i = 1, n_index_nz(nz)
        iii = iii + 1

        write (2,*) nodest(1,iii),nodest(2,iii),nodest(3,iii),nodest(4,iii)
        write (2,*) hostst(1,iii),hostst(2,iii),hostst(3,iii),hostst(4,iii)
        write (2,*)  coeft(1,iii), coeft(2,iii), coeft(3,iii)
       if(coeft(1,i)>1. .or. coeft(2,i)>1. .or. coeft(3,i)>1.) then
          print *, 'coef', i,nz, nodest(1:4,i)
          print *, 'node', i,nz, coeft(1:4,i)
          end if
     end do
     print *, 'sum(coeft)',nz,sum(coeft(:,:))
  end do
  close(2)


  
END PROGRAM MAIN

SUBROUTINE search(point,cell,inside,dis)
  implicit none
  real::point(3),cell(8,3),dis(6)
  integer::inside,i,p1(6),p2(6),p3(6),p4(6)
  inside=0

  p1(5)=1
  p2(5)=4
  p3(5)=3
  p4(5)=2

  p1(3)=1
  p2(3)=2
  p3(3)=6
  p4(3)=5

  p1(1)=1
  p2(1)=5
  p3(1)=8
  p4(1)=4

  p1(2)=2
  p2(2)=3
  p3(2)=7
  p4(2)=6

  p1(6)=5
  p2(6)=6
  p3(6)=7
  p4(6)=8

  p1(4)=4
  p2(4)=8
  p3(4)=7
  p4(4)=3

  do i=1,6
     call vector(cell(p1(i),:),cell(p2(i),:),cell(p3(i),:),cell(p4(i),:),  &
          &	point,dis(i))
     !if(abs(dis(i))<1.0E-06) dis(i)=0.

     if(dis(i).lt.0.0)then
        inside=0
        exit
     else
        inside=1
     end if

  end do

END SUBROUTINE search

SUBROUTINE vector(p1,p2,p3,p4,point,dis)
  implicit none
  real::dis,p1(3),p2(3),p3(3),p4(3),point(3)
  real::a1v,a2v,a3v,b1v,b2v,b3v,c1v,c2v,c3v,ab1,ab2,ab3

  a1v=p4(1)-p2(1)
  a2v=p4(2)-p2(2)
  a3v=p4(3)-p2(3)
  
  b1v=p3(1)-p1(1)
  b2v=p3(2)-p1(2)
  b3v=p3(3)-p1(3)

  c1v=point(1)-(p1(1)+p2(1)+p3(1)+p4(1))/4.
  c2v=point(2)-(p1(2)+p2(2)+p3(2)+p4(2))/4.
  c3v=point(3)-(p1(3)+p2(3)+p3(3)+p4(3))/4.

  ab1=a2v*b3v-a3v*b2v
  ab2=a3v*b1v-a1v*b3v
  ab3=a1v*b2v-a2v*b1v

  dis=c1v*ab1+c2v*ab2+c3v*ab3
END SUBROUTINE vector

!c--------------------------------------------------------------------
!c
!c       A globally convergent Newton-Raphson method (from Numerical
!c       receipes in Fortran, by  W. H. Press, 2nd ed., p379 ) 
!c       Parameter: NP--dimension of the sysytem, MAXITS--maximum 
!c       number of iterations within which convergence shall be 
!c       obtained if it is supposed. 
!c
!c-------------------------------------------------------------------


SUBROUTINE newt(x,n,check,MAXITS,its,itss)
  
  INTEGER n,nn,NP
  LOGICAL check
  REAL x(n),fvec,TOLF,TOLMIN,TOLX,STPMX

  PARAMETER (NP=40,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-7,	&
       &  STPMX=100.)
  COMMON /newtv/ fvec(NP),nn
  SAVE /newtv/
!C     USES fdjac,fmin,lnsrch,lubksb,ludcmp

  INTEGER i,j,indx(NP),MAXITS
  REAL d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),p(NP),	&
       &  xold(NP),fmin
  EXTERNAL fmin
  nn=n
  f=fmin(x)
  test=0.
  do i=1,n
     if(abs(fvec(i)).gt.test)test=abs(fvec(i))
  end do
  if(test.lt..01*TOLF) then
     check=.false.
     return
  endif
  sum=0.
  do i=1,n
     sum=sum+x(i)**2
  end do
  stpmax=STPMX*max(sqrt(sum),float(n))
  do its=1,MAXITS
     write(*,*) its
!     call fdjac(n,x,fvec,NP,fjac)
     call fdjac1(n,x,fjac(1:n,1:n))
     call funcv(n,x,fvec)
!     print *,fjac(1:n,1:n)
     do i=1,n
        sum=0.
        do j=1,n
           sum=sum+fjac(j,i)*fvec(j)
        end do
        g(i)=sum
     end do
     do i=1,n
        xold(i)=x(i)
     end do
     fold=f
     do i=1,n
        p(i)=-fvec(i)
     end do
!     write(*,*) '********'
     call ludcmp(fjac,n,NP,indx,d,itss)
!     write(*,*) '*******1'
     call lubksb(fjac,n,NP,indx,p)
!     write(*,*) '*******2'
     call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
!     write(*,*) '*******3'
     test=0.
     do i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
     end do
     if(test.lt.TOLF)then
        check=.false.
        return
     endif
     if(check)then
        test=0.
        den=max(f,.5*n)
        do i=1,n
           temp=abs(g(i))*max(abs(x(i)),1.)/den
           if(temp.gt.test)test=temp
        end do
        if(test.lt.TOLMIN)then
           check=.true.
        else
           check=.false.
        endif
        return
     endif
     test=0.
     do i=1,n
        temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
        if(temp.gt.test)test=temp
     end do
     if(test.lt.TOLX)return
  end do
!c      pause 'MAXITS exceeded in newt'
END SUBROUTINE newt


SUBROUTINE fdjac(n,x,fvec,np,df)
  INTEGER n,np,NMAX
  REAL df(np,np),fvec(n),x(n),EPS
  PARAMETER (NMAX=40,EPS=1.e-4)
!    USES funcv
  INTEGER i,j
  REAL h,temp,f(NMAX)
  do j=1,n
     temp=x(j)
     h=EPS*abs(temp)
     if(h.eq.0.)h=EPS
     x(j)=temp+h
     h=x(j)-temp
     call funcv(n,x,f)
     x(j)=temp
     do i=1,n
        df(i,j)=(f(i)-fvec(i))/h
     end do
  end do
  return
END SUBROUTINE fdjac


FUNCTION fmin(x)
  INTEGER n,NP
  REAL fmin,x(*),fvec
  PARAMETER (NP=40)
  COMMON /newtv/ fvec(NP),n
  SAVE /newtv/
!CU    USES funcv
  INTEGER i
  REAL sum
  call funcv(n,x,fvec)
  sum=0.
  do i=1,n
     sum=sum+fvec(i)**2
  end do
  fmin=0.5*sum
  return
END FUNCTION fmin

SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
  INTEGER n
  LOGICAL check
  REAL f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
  PARAMETER (ALF=1.e-4,TOLX=1.e-7)
  EXTERNAL func
!CU    USES func
  INTEGER i
  REAL a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,	&
       &  test,tmplam
  check=.false.
  sum=0.
  do i=1,n
     sum=sum+p(i)*p(i)
  end do
  sum=sqrt(sum)
  if(sum.gt.stpmax)then
     do i=1,n
        p(i)=p(i)*stpmax/sum
     end do
  endif
  slope=0.
  do i=1,n
     slope=slope+g(i)*p(i)
  end do
  if (slope.ge.0.) pause 'roundoff problem in lnsrch'
  test=0.
  do i=1,n
     temp=abs(p(i))/max(abs(xold(i)),1.)
     if(temp.gt.test)test=temp
  end do
  alamin=TOLX/test
  alam=1.

  alam2=2.
  f2=2.
  fold2=fold

1     continue
  do i=1,n
     x(i)=xold(i)+alam*p(i)
  end do
  f=func(x)
  if(alam.lt.alamin)then
     do i=1,n
        x(i)=xold(i)
     end do
     check=.true.
     return
  else if(f.le.fold+ALF*alam*slope)then
     return
  else
     if(alam.eq.1.)then
        tmplam=-slope/(2.*(f-fold-slope))
     else
        rhs1=f-fold-alam*slope
        rhs2=f2-fold-alam2*slope
        a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
        b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
        if(a.eq.0.)then
           tmplam=-slope/(2.*b)
        else
           disc=b*b-3.*a*slope
           if(disc.lt.0.) then
              tmplam=0.5*alam
           else if (b.le.0.) then
              tmplam=(-b+sqrt(disc))/(3.*a)
           else
              tmplam=-slope/(b+sqrt(disc))
           endif
!           C              tmplam=(-b+sqrt(disc))/(3.*a)
        endif
        if(tmplam.gt..5*alam)tmplam=.5*alam
     endif
  endif
  alam2=alam
  f2=f
!  fold2=fold
  alam=max(tmplam,.1*alam)
  goto 1
END SUBROUTINE lnsrch

SUBROUTINE ludcmp(a,n,np,indx,d,itss)
  INTEGER n,np,indx(n),NMAX
  REAL d,a(np,np),TINY
  PARAMETER (NMAX=500,TINY=1.0e-20)
  INTEGER i,imax,j,k
  REAL aamax,dum,sum,vv(NMAX)
  itss=0
  d=1.
  do i=1,n
     aamax=0.
     do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
     end do
     if (aamax.eq.0.) pause 'singular matrix in ludcmp'
     if (aamax.eq.0.) itss=1
     if (aamax.eq.0.) return
     vv(i)=1./aamax
  end do
  do j=1,n
     do i=1,j-1
        sum=a(i,j)
        do k=1,i-1
           sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
     end do
     aamax=0.
     do i=j,n
        sum=a(i,j)
        do k=1,j-1
           sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
        dum=vv(i)*abs(sum)
        if (dum.ge.aamax) then
           imax=i
           aamax=dum
        endif
     end do
     if (j.ne.imax)then
        do k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        end do
        d=-d
        vv(imax)=vv(j)
     endif
     indx(j)=imax
     if(a(j,j).eq.0.)a(j,j)=TINY
     if(j.ne.n)then
        dum=1./a(j,j)
        do i=j+1,n
           a(i,j)=a(i,j)*dum
        end do
     endif
  end do
  return
END SUBROUTINE ludcmp

SUBROUTINE lubksb(a,n,np,indx,b)
  INTEGER n,np,indx(n)
  REAL a(np,np),b(n)
  INTEGER i,ii,j,ll
  REAL sum
  ii=0
  do i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii.ne.0)then
        do j=ii,i-1
           sum=sum-a(i,j)*b(j)
        end do
     else if (sum.ne.0.) then
        ii=i
     endif
     b(i)=sum
  end do
  do i=n,1,-1
     sum=b(i)
     do j=i+1,n
        sum=sum-a(i,j)*b(j)
     end do
     b(i)=sum/a(i,i)
  end do
  return
END SUBROUTINE lubksb

SUBROUTINE funcv(n,rst,fvec)

  real::p(3),v000(3),v100(3),v010(3),v001(3),v101(3),v011(3),v110(3),v111(3)
  common/intp/v000,v100,v010,v001,v101,v011,	&
       &   v110,v111,p
  real::rst(n),fvec(n)
  do i=1,n
     fvec(i)=v000(i)*(1-rst(1))*(1-rst(2))*(1-rst(3))+	&
          &  v100(i)*rst(1)*(1-rst(2))*(1-rst(3))+	&
          &  v010(i)*(1-rst(1))*rst(2)*(1-rst(3))+	&
          &  v001(i)*(1-rst(1))*(1-rst(2))*rst(3)+	&
          &  v101(i)*rst(1)*(1-rst(2))*rst(3)+		&
          &  v011(i)*(1-rst(1))*rst(2)*rst(3)+		&
          &  v110(i)*rst(1)*rst(2)*(1-rst(3))+		&
          &  v111(i)*rst(1)*rst(2)*rst(3)-p(i)
  end do
  return
END SUBROUTINE funcv

SUBROUTINE fdjac1(n,rst,fjac)
  real::p(3),v000(3),v100(3),v010(3),v001(3),v101(3),v011(3),v110(3),v111(3)
  common/intp/v000,v100,v010,v001,v101,v011,	&
       &   v110,v111,p
  real::rst(n),fjac(n,n)

  fjac(1,1)=-v000(1)*(1-rst(2))*(1-rst(3))+	&
       &  v100(1)*(1-rst(2))*(1-rst(3))-	&
       &  v010(1)*rst(2)*(1-rst(3))-	&
       &  v001(1)*(1-rst(2))*rst(3)+	&
       &  v101(1)*(1-rst(2))*rst(3)-	&
       &  v011(1)*rst(2)*rst(3)+	&
       &  v110(1)*rst(2)*(1-rst(3))+	&
       &  v111(1)*rst(2)*rst(3)
  fjac(2,1)=-v000(2)*(1-rst(2))*(1-rst(3))+	&
       &  v100(2)*(1-rst(2))*(1-rst(3))-	&
       &  v010(2)*rst(2)*(1-rst(3))-	&
       &  v001(2)*(1-rst(2))*rst(3)+	&
       &  v101(2)*(1-rst(2))*rst(3)-	&
       &  v011(2)*rst(2)*rst(3)+	&
       &  v110(2)*rst(2)*(1-rst(3))+	&
       &  v111(2)*rst(2)*rst(3)
  fjac(3,1)=-v000(3)*(1-rst(2))*(1-rst(3))+	&
       &  v100(3)*(1-rst(2))*(1-rst(3))-	&
       &  v010(3)*rst(2)*(1-rst(3))-	&
       &  v001(3)*(1-rst(2))*rst(3)+	&
       &  v101(3)*(1-rst(2))*rst(3)-	&
       &  v011(3)*rst(2)*rst(3)+	&
       &  v110(3)*rst(2)*(1-rst(3))+	&
       &  v111(3)*rst(2)*rst(3)
  fjac(1,2)=-v000(1)*(1-rst(1))*(1-rst(3))-	&
       &  v100(1)*rst(1)*(1-rst(3))+	&
       &  v010(1)*(1-rst(1))*(1-rst(3))-	&
       &  v001(1)*(1-rst(1))*rst(3)-	&
       &  v101(1)*rst(1)*rst(3)+	&
       &  v011(1)*(1-rst(1))*rst(3)+	&
       &  v110(1)*rst(1)*(1-rst(3))+	&
       &  v111(1)*rst(1)*rst(3)
  fjac(2,2)=-v000(2)*(1-rst(1))*(1-rst(3))-	&
       &  v100(2)*rst(1)*(1-rst(3))+	&
       &  v010(2)*(1-rst(1))*(1-rst(3))-	&
       &  v001(2)*(1-rst(1))*rst(3)-	&
       &  v101(2)*rst(1)*rst(3)+	&
       &  v011(2)*(1-rst(1))*rst(3)+	&
       &  v110(2)*rst(1)*(1-rst(3))+	&
       &  v111(2)*rst(1)*rst(3)
  fjac(3,2)=-v000(3)*(1-rst(1))*(1-rst(3))-	&
       &  v100(3)*rst(1)*(1-rst(3))+	&
       &  v010(3)*(1-rst(1))*(1-rst(3))-	&
       &  v001(3)*(1-rst(1))*rst(3)-	&
       &  v101(3)*rst(1)*rst(3)+	&
       &  v011(3)*(1-rst(1))*rst(3)+	&
       &  v110(3)*rst(1)*(1-rst(3))+	&
       &  v111(3)*rst(1)*rst(3)
  fjac(1,3)=-v000(1)*(1-rst(1))*(1-rst(2))-	&
       &  v100(1)*rst(1)*(1-rst(2))-	&
       &  v010(1)*(1-rst(1))*rst(2)+	&
       &  v001(1)*(1-rst(1))*(1-rst(2))+	&
       &  v101(1)*rst(1)*(1-rst(2))+	&
       &  v011(1)*(1-rst(1))*rst(2)-	&
       &  v110(1)*rst(1)*rst(2)+	&
       &  v111(1)*rst(1)*rst(2)
  fjac(2,3)=-v000(2)*(1-rst(1))*(1-rst(2))-	&
       &  v100(2)*rst(1)*(1-rst(2))-	&
       &  v010(2)*(1-rst(1))*rst(2)+	&
       &  v001(2)*(1-rst(1))*(1-rst(2))+	&
       &  v101(2)*rst(1)*(1-rst(2))+	&
       &  v011(2)*(1-rst(1))*rst(2)-	&
       &  v110(2)*rst(1)*rst(2)+	&
       &  v111(2)*rst(1)*rst(2)
  fjac(3,3)=-v000(3)*(1-rst(1))*(1-rst(2))-	&
       &  v100(3)*rst(1)*(1-rst(2))-	&
       &  v010(3)*(1-rst(1))*rst(2)+	&
       &  v001(3)*(1-rst(1))*(1-rst(2))+	&
       &  v101(3)*rst(1)*(1-rst(2))+	&
       &  v011(3)*(1-rst(1))*rst(2)-	&
       &  v110(3)*rst(1)*rst(2)+	&
       &  v111(3)*rst(1)*rst(2)
END SUBROUTINE fdjac1
