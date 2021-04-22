  MODULE  emparam
      use prec
      implicit none
      real(dp) :: delta
      character*1 :: icf
      real(dp), dimension(2,3) :: kvcbm
      integer :: nkvc, np
  END MODULE emparam

  subroutine readinpem()
   use prec 
   use emparam

   implicit none
   integer :: i, j , k
   open(40,file='inp.em',status='old')
   read(40,*) (kvcbm(1,i),i=1,3), (kvcbm(2,i),i=1,3)
   read(40,*) icf
   read(40,*) delta , np
   close(40)

   if ( (abs(kvcbm(1,1) - kvcbm(2,1)) .LT. 0.001d0) .and.  &
&       (abs(kvcbm(1,2) - kvcbm(2,2)) .LT. 0.001d0) .and.  &
&       (abs(kvcbm(1,3) - kvcbm(2,3)) .LT. 0.001d0) ) then
     nkvc = 1
   else 
     nkvc = 2
   end if
  
   end subroutine readinpem

    subroutine gen_kline_around_kp(kp, id, n, delta, kptlines)
    use prec
    implicit none
    integer:: i, j
    integer,  intent(in), dimension(3)::   id
    integer,  intent(in)::   n
    real(dp), intent(in), dimension(3)::   kp
    real(dp), intent(in)::   delta
    real(dp), intent(out), dimension(2*n+1,3)::  kptlines
    j=0
    do  i= -n, n
       if ( i .NE. 0) then
          j= j+1
          kptlines(j,1) =  kp(1) + i * dble(id(1)) * delta
          kptlines(j,2) =  kp(2) + i * dble(id(2)) * delta
          kptlines(j,3) =  kp(3) + i * dble(id(3)) * delta
      endif
    end do
    return
    end subroutine gen_kline_around_kp

    subroutine gen_kgrid_around_kp(kp, n, delta, kptg)
    use prec
    implicit none
    integer:: i, j, k, ik
    integer, intent(in)::   n
    real(dp), intent(in), dimension(3) ::   kp
    real(dp), intent(in) :: delta
    real(dp), intent(out), dimension((2*n+1)**3,3):: kptg
    ik=0
    do i = -n, n
        do j = -n, n
          do k = -n, n
            ik= ik+1
            kptg(ik,1) = kp(1)+ dble(i) * delta 
            kptg(ik,2) = kp(2)+ dble(j) * delta 
            kptg(ik,3) = kp(3)+ dble(k) * delta 
          end do  
        end do
    end do 

    return
    end subroutine gen_kgrid_around_kp


