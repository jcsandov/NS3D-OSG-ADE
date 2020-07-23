
subroutine tracking ()

  implicit none

  ! local variables
integer ::nn
integer,dimension(5)::n_int
integer,dimension(5,4)::zint

n_int(1) = 0 !2
n_int(2) = 1 !2
n_int(3) = 0 !3
n_int(4) = 1 !3
n_int(5) = 0 !3

zint(1,1) = 5; zint(1,2) = 4
zint(2,1) = 4; zint(2,2) = 4
zint(3,1) = 2; zint(3,2) = 2; zint(3,3) = 4
zint(4,1) = 1; zint(4,2) = 3; zint(4,3) = 5
zint(5,1) = 1; zint(5,2) = 3; zint(5,3) = 4




! START TRACKING
!
nz = szone(n)


ibck = ii(n)-3
ifwd = ii(n)+3
jbck = jj(n)-3
jfwd = jj(n)+3
kbck = kk(n)-3
kfwd = kk(n)+3

iflag = 0

call locate(nz, im(nz), jm(nz), km(nz), &
	x_l(1:im(nz),1:jm(nz),1:km(nz),nz), &
	y_l(1:im(nz),1:jm(nz),1:km(nz),nz), &
	z_l(1:im(nz),1:jm(nz),1:km(nz),nz), &
	qtotal(1:3,1:im(nz),1:jm(nz),1:km(nz),nz), &
	xp(n), yp(n), zp(n), ii(n), jj(n), kk(n), &
	ibck, ifwd, jbck, jfwd, kbck, kfwd, &
	u_vel(n), v_vel(n), w_vel(n), &
	Du(1:3,1:im(nz),1:jm(nz),1:km(nz),nz), &
	vort(1:3,1:im(nz),1:jm(nz),1:km(nz),nz), &
       Dux(n), Duy(n), Duz(n), &
      vortx(n), vorty(n), vortz(n), &
	iflag, status_n(n),treeflag)



if (iflag.eq.0) then
nz = szone(n)
treeflag = 1
ibck = 1; ifwd = im(nz) - 1
jbck = 1; jfwd = jm(nz) - 1
kbck = 1; kfwd = km(nz) - 1

	call locate(nz, im(nz), jm(nz), km(nz), &
		x_l(1:im(nz),1:jm(nz),1:km(nz),nz), &
		y_l(1:im(nz),1:jm(nz),1:km(nz),nz), &
		z_l(1:im(nz),1:jm(nz),1:km(nz),nz), &
		qtotal(1:3,1:im(nz),1:jm(nz),1:km(nz),nz), &
		xp(n), yp(n), zp(n), ii(n), jj(n), kk(n), &
		ibck, ifwd, jbck, jfwd, kbck, kfwd, &
		u_vel(n), v_vel(n), w_vel(n), &
		Du(1:3,1:im(nz),1:jm(nz),1:km(nz),nz), &
		vort(1:3,1:im(nz),1:jm(nz),1:km(nz),nz), &
     	  	Dux(n), Duy(n), Duz(n), &
     	 	vortx(n), vorty(n), vortz(n), &
		iflag, status_n(n),treeflag)

end if



! if particle is not found, change zone
if (iflag.eq.0) then
 	treeflag = 1

	do nn=1,n_int(szone(n))
		nz=zint(szone(n),nn)

		ibck = 1
		ifwd = im(nz) - 1
		jbck = 1
		jfwd = jm(nz) - 1
		kbck = 1
		kfwd = km(nz) - 1


	call locate(nz, im(nz), jm(nz), km(nz), &
		x_l(1:im(nz),1:jm(nz),1:km(nz),nz), &
		y_l(1:im(nz),1:jm(nz),1:km(nz),nz), &
		z_l(1:im(nz),1:jm(nz),1:km(nz),nz), &
		qtotal(1:3,1:im(nz),1:jm(nz),1:km(nz),nz), &
		xp(n), yp(n), zp(n), ii(n), jj(n), kk(n), &
		ibck, ifwd, jbck, jfwd, kbck, kfwd, &
		u_vel(n), v_vel(n), w_vel(n), &
		Du(1:3,1:im(nz),1:jm(nz),1:km(nz),nz), &
		vort(1:3,1:im(nz),1:jm(nz),1:km(nz),nz), &
     	 	Dux(n), Duy(n), Duz(n), &
     	 	vortx(n), vorty(n), vortz(n), &
		iflag, status_n(n), treeflag )

		if (iflag.eq.1) then
			szone(n) = nz
			treeflag = 0
			return
		end if

	end do
end if



! if particle is not found, it exited the domain
! and it won't be calculated any longer.
if (iflag.eq.0) then
	status_n(n) = 0
      print *,'particle=',n,'in processor=',myid,'exited'
      print *,'x,y,z,szone',xp(n),yp(n),zp(n),szone(n)
      print *,'last known position',ii(n),jj(n),kk(n)
end if

treeflag = 0

end subroutine tracking

