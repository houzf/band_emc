   SUBROUTINE get_kdistances(n,kpts,icf,kdist)
     use prec
     use latticevector

     implicit none
     integer i, j
     integer, intent(in):: n
     character(len=1), intent(in) ::  icf
     real(dp), intent(in), dimension(n,3) :: kpts 
     real(dp), intent(out), dimension(n) :: kdist 
     real(dp), dimension(3) :: kbuf

     call readposcar()
     kdist(1) = 0.0d0

     if ((icf .EQ. 'C') .OR. (icf .EQ. 'c')) then
       do i = 2, n
          kbuf=kpts(i,:) - kpts(i-1,:) 
          kdist(i)= kdist(i-1) + sqrt(DOT_PRODUCT(kbuf, kbuf)) 
       end do
     else
       do i = 2, n
          kbuf=MATMUL(kpts(i,:),recip_lat) - &
 &               MATMUL(kpts(i-1,:),recip_lat) 
          kdist(i)= kdist(i-1) + sqrt(DOT_PRODUCT(kbuf, kbuf)) 
       end do
     end if

    end subroutine get_kdistances 
