
subroutine rhs_ns_unst_visc_diss


  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Calculate second-order accurate unsteady terms and add
  ! current viscous and dissipation terms to rh

  ! Unlike previous codes we will add viscous, dissipation, and
  ! unsteady terms in the calling routine not in this one.

  ! input
  !     aj(ijk)
  !     q(4,ijk)
  !     qn(4,ijk)
  !     qnm1(4,ijk)
  !     delti
  !     rh(4,ijk)

  ! output
  !     rh(4,ijk)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! local dummy reals
  real (kind = rdf) :: dtaj
  real (kind = rdf) :: dq2dt, dq3dt, dq4dt

  ! Add time derivative, viscous and dissipation terms to rh

   if (nlinc) then

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend

     dtaj = two * delti * aj(i,j,k)
     dq2dt = e_source*(three*q(2,i,j,k)-four*qn(2,i,j,k)+qnm1(2,i,j,k)) / dtaj
     dq3dt = e_source*(three*q(3,i,j,k)-four*qn(3,i,j,k)+qnm1(3,i,j,k)) / dtaj
     dq4dt = e_source*(three*q(4,i,j,k)-four*qn(4,i,j,k)+qnm1(4,i,j,k)) / dtaj

     rh(1,i,j,k) = rh(1,i,j,k) + diss(1,i,j,k)
     rh(2,i,j,k) = rh(2,i,j,k) + diss(2,i,j,k) - visc(1,i,j,k) + dq2dt + &
                   rstr(1,i,j,k)
     rh(3,i,j,k) = rh(3,i,j,k) + diss(3,i,j,k) - visc(2,i,j,k) + dq3dt + &
                   rstr(2,i,j,k)
     rh(4,i,j,k) = rh(4,i,j,k) + diss(4,i,j,k) - visc(3,i,j,k) + dq4dt + &
                   rstr(3,i,j,k)

   end do
   end do
   end do

   else

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend

     dtaj = two * delti * aj(i,j,k)
     dq2dt = e_source*(three*q(2,i,j,k)-four*qn(2,i,j,k)+qnm1(2,i,j,k)) / dtaj
     dq3dt = e_source*(three*q(3,i,j,k)-four*qn(3,i,j,k)+qnm1(3,i,j,k)) / dtaj
     dq4dt = e_source*(three*q(4,i,j,k)-four*qn(4,i,j,k)+qnm1(4,i,j,k)) / dtaj

     rh(1,i,j,k) = rh(1,i,j,k) + diss(1,i,j,k)
     rh(2,i,j,k) = rh(2,i,j,k) + diss(2,i,j,k) - visc(1,i,j,k) + dq2dt
     rh(3,i,j,k) = rh(3,i,j,k) + diss(3,i,j,k) - visc(2,i,j,k) + dq3dt
     rh(4,i,j,k) = rh(4,i,j,k) + diss(4,i,j,k) - visc(3,i,j,k) + dq4dt

   end do
   end do
   end do

   end if

end subroutine rhs_ns_unst_visc_diss

