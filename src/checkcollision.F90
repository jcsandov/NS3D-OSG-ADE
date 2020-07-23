
subroutine checkcollision ()

! this collision algorithm is still
! too simplified.  It works as long as the
! time step is small. A collision detection
! scheme will be required in the future.

  implicit none

  ! local variables
!integer ::nn
real (kind = rdf):: cost, sint, radip


if ( zp(n).lt.(ds(n)/two).and.wp(n).lt.zero) then

zp(n) = ds(n) / two 
up(n) = 0.9*up(n)
vp(n) = 0.9*vp(n)
wp(n) = zero !-wp(n)

                 xold(1) = xp(n)
                 xold(2) = yp(n)
                 xold(3) = zp(n)

                 xold(4) = up(n)
                 xold(5) = vp(n)
                 xold(6) = wp(n)

hitflag = 1

end if


if ( zp(n).lt.zero) then

zp(n) = ds(n) / two 
up(n) = 0.9*up(n)
vp(n) = 0.9*vp(n)
wp(n) = zero !-wp(n)

                 xold(1) = xp(n)
                 xold(2) = yp(n)
                 xold(3) = zp(n)

                 xold(4) = up(n)
                 xold(5) = vp(n)
                 xold(6) = wp(n)

hitflag = 1

end if



if (sqrt(xp(n)*xp(n)+(yp(n)-five)*(yp(n)-five)).lt.pt5+ds(n)/two) then

	radip = pt5 + ds(n)/two
	xp(n) = sign(sqrt(radip**two - (yp(n)-five)**two), xp(n))
	cost = xp(n)/radip
	sint = (yp(n)-five)/radip
	up(n) =  0.9*sint*(up(n)*sint-vp(n)*cost)
	vp(n) = -0.9*cost*(up(n)*sint-vp(n)*cost)
	wp(n) =  0.9*wp(n)
	
                 xold(1) = xp(n)
                 xold(2) = yp(n)
                 xold(3) = zp(n)

                 xold(4) = up(n)
                 xold(5) = vp(n)
                 xold(6) = wp(n)

	hitflag = 1

end if


end subroutine checkcollision

