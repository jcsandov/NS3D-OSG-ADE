
subroutine bcond_fm (il, iu, jl, ju, kl, ku, igp, jgp, kgp, q)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  ! general version !Kim's curved duct (full 3d geometry)
  ! 
  ! fine-mesh
  ! boundary conditions
  ! 
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  use global_app
  use global_mpi

  implicit none
  
  ! extents
  ! 
  integer :: il
  integer :: jl
  integer :: kl
  integer :: iu
  integer :: ju
  integer :: ku

  ! ghost points
  ! 
  integer :: igp
  integer :: jgp
  integer :: kgp

  ! solution vector
  ! 
  real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku) :: q

  ! fixed pressure
  ! 
  real (kind = rdf) :: pfix

  integer :: i, j, k, b
  
  integer :: i_mysta
  integer :: j_mysta
  integer :: k_mysta

  integer :: i_myend
  integer :: j_myend
  integer :: k_myend

  ! make definition of i_mysta etc., consistent with
  ! all other uses, e.g., solver_daf & rhs_kw_bcond
  !
  i_mysta = il + igp
  j_mysta = jl + jgp
  k_mysta = kl + kgp

  i_myend = iu - igp
  j_myend = ju - jgp
  k_myend = ku - kgp

  ! processes on the domain boundaries
  ! 
  !if (myback == mpi_proc_null)  i_mysta = il + igp + 1
  !if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
  !if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

  !if (myfront == mpi_proc_null) i_myend = iu - igp - 1
  !if (myright == mpi_proc_null) j_myend = ju - jgp - 1
  !if (myup    == mpi_proc_null) k_myend = ku - kgp - 1


  ! boundary in csi-direction (i = 1)
  !
  if (myback == mpi_proc_null) then
     b = 1
     i = i_mysta

     csi1: select case(btype(b,myzone))
     case(0) ! -> interface

     case(1:3) ! wall, symmetric plane & freestream
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i+1,j,k) + &
                          sb(1:4,b) * q(1:4,i+2,j,k)
        end do
        end do
     case(4) ! inflow
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           q(1,i,j,k) = (one + rat(j,k,b)) * q(1,i+1,j,k) - &
                               rat(j,k,b)  * q(1,i+2,j,k)
           !q(1,i,j,k) = sa(1,b) * q(1,i+1,j,k) + &
           !             sb(1,b) * q(1,i+2,j,k)
        end do
        end do
     case(5) ! exit (MOC)
     case(6) ! periodic condition
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           q(1:4,i,j,k) = q(1:4,i_myend,j,k)
        end do
        end do
     end select csi1
  end if
     
  ! boundary in csi-direction (i = imax)
  !
  if (myfront == mpi_proc_null) then
     b = 2
     i = i_myend

     csi2: select case(btype(b,myzone))
     case(0) ! -> interface
     case(1:3) ! wall, symmetric plane & freestream
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i-1,j,k) + &
                          sb(1:4,b) * q(1:4,i-2,j,k)
        end do
        end do
     case(4) ! inflow
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           q(1,i,j,k) = (one + rat(j,k,b)) * q(1,i-1,j,k) - &
                               rat(j,k,b)  * q(1,i-2,j,k)
           !q(1,i,j,k) = sa(1,b) * q(1,i-1,j,k) + &
           !             sb(1,b) * q(1,i-2,j,k)
        end do
        end do
     case(5) ! exit (MOC)
		!j = j_mysta
        !do k = k_mysta, k_myend
		!   q( 1 ,i,j,k) = q( 1 ,i-1,j,k)
        !   q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i-1,j,k) + &
        !                  sb(2:4,b) * q(2:4,i-2,j,k)
		!end do

		!j = j_myend
        !do k = k_mysta, k_myend
		!   q( 1 ,i,j,k) = q( 1 ,i-1,j,k)
        !   q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i-1,j,k) + &
        !                  sb(2:4,b) * q(2:4,i-2,j,k)
		!end do

     case(6) ! periodic condition
        !do k = k_mysta, k_myend
        !do j = j_mysta, j_myend
        !   q(1:4,i,j,k) = q(1:4,i_mysta,j,k)
        !end do
        !end do
     case(8) ! free surface
          do k = k_mysta, k_myend
          do j = j_mysta, j_myend
           q(1:4,i,j,k) = q(1:4,i-1,j,k) !sa(1:4,b)*q(1:4,i-1,j,k) +
                                         !sb(1:4,b)*q(1:4,i-2,j,k)                
          end do
          end do

    end select csi2
  end if
     
  ! boundary in eta-direction (j = 1)
  ! 
  if (myleft == mpi_proc_null) then
     b = 3
     j = j_mysta

     eta1: select case(btype(b,myzone))
     case(0) ! -> interface
     case(1:3) ! wall, symmetric plane & freestream
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j+1,k) + &
                          sb(1:4,b) * q(1:4,i,j+2,k)
        end do
        end do
     case(4) ! inflow
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1,i,j,k) = (one + rat(i,k,b)) * q(1,i,j+1,k) - &
                               rat(i,k,b)  * q(1,i,j+2,k)
           !q(1,i,j,k) = sa(1,b) * q(1,i,j+1,k) + &
           !             sb(1,b) * q(1,i,j+2,k)
        end do
        end do
     case(5) ! exit (MOC)
     case(6) ! periodic condition
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = q(1:4,i,j_myend,k)
        end do
        end do
     end select eta1
  end if

  ! boundary in eta-direction (j = jmax)
  ! 
  if (myright == mpi_proc_null) then
     b = 4
     j = j_myend

     eta2: select case(btype(b,myzone))
     case(0) ! -> interface
     case(1:3) ! wall, symmetric plane & freestream
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j-1,k) + &
                          sb(1:4,b) * q(1:4,i,j-2,k)
        end do
        end do
     case(4) ! inflow
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1,i,j,k) = (one + rat(i,k,b)) * q(1,i,j-1,k) - &
                               rat(i,k,b)  * q(1,i,j-2,k)
           !q(1,i,j,k) = sa(1,b) * q(1,i,j-1,k) + &
           !             sb(1,b) * q(1,i,j-2,k)
        end do
        end do
     case(5) ! exit (MOC)
     case(6) ! periodic condition
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = q(1:4,i,j_mysta,k)
        end do
        end do
     end select eta2
  end if

  ! boundary in zet-direction (k = 1)
  ! 
  if (mydown == mpi_proc_null) then
     b = 5
     k = k_mysta

     zet1: select case(btype(b,myzone))
     case(0) ! -> interface
     case(1:3) ! wall, symmetric plane & freestream
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j,k+1) + &
                          sb(1:4,b) * q(1:4,i,j,k+2)
        end do
        end do
     case(4) ! inflow
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1,i,j,k) = (one + rat(i,j,b)) * q(1,i,j,k+1) - &
                               rat(i,j,b)  * q(1,i,j,k+2)
           !q(1,i,j,k) = sa(1,b) * q(1,i,j,k+1) + &
           !             sb(1,b) * q(1,i,j,k+2)
        end do
        end do
     case(5) ! exit (MOC)
     case(6) ! periodic condition
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = q(1:4,i,j,k_myend)
        end do
        end do
     end select zet1
  end if

  ! boundary in zet-direction (k = kmax)
  ! 
  if (myup == mpi_proc_null) then
     b = 6
     k = k_myend

     zet2: select case(btype(b,myzone))
     case(0) ! -> interface
     case(1:3) ! wall, symmetric plane & freestream
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j,k-1) + &
                          sb(1:4,b) * q(1:4,i,j,k-2)
        end do
        end do
     case(4) ! inflow
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1,i,j,k) = (one + rat(i,j,b)) * q(1,i,j,k-1) - &
                               rat(i,j,b)  * q(1,i,j,k-2)
           !q(1,i,j,k) = sa(1,b) * q(1,i,j,k-1) + &
           !             sb(1,b) * q(1,i,j,k-2)
        end do
        end do
     case(5) ! exit (MOC)
     case(6) ! periodic condition
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = q(1:4,i,j,k_mysta)
        end do
        end do
     end select zet2
  end if

  ! probe default pressure value
  !
  if (myid == pfix_proc) pfix = q(1,local_ifix,local_jfix,local_kfix)

  ! distribute fixed pressure to all processes
  ! 
  call mpi_bcast(pfix, 1, mpi_real, pfix_proc, mpi_comm_world, ierr)

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     q(1,i,j,k) = q(1,i,j,k) - pfix
  end do
  end do
  end do

end subroutine bcond_fm

