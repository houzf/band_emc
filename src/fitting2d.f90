! Adapted for 2-dimensional and 3-dimensional polynomial fitting
! http://rosettacode.org/wiki/Polynomial_Fitting#Fortran
module fitting2d
contains
 
  function polyfit(vx, vy, d, dimens)
    implicit none
    integer, intent(in)                   :: d, dimens
    integer, parameter                    :: dp = selected_real_kind(15, 307)
    real(dp), dimension((d+1)*(d+2)*(d+3)/6)  :: polyfit
    real(dp), dimension(:,:), intent(in)  :: vx
    real(dp), dimension(:), intent(in)    :: vy
 
    real(dp), dimension(:,:), allocatable :: X
    real(dp), dimension(:,:), allocatable :: XT
    real(dp), dimension(:,:), allocatable :: XTX
 
    integer :: i, j, k, ip, iq
 
    integer     :: n, lda, lwork
    integer :: info
    integer, dimension(:), allocatable :: ipiv
    real(dp), dimension(:), allocatable :: work
 
    SELECT CASE (dimens)
        case (1)
          n = d + 1
        case (2)
          n = (d + 1) * (d + 2) / 2
        case (3)
          n = (d + 1) * (d + 2) * (d + 3) / 6
        case default
             print *,  "Not implemented for 4 or more variables"
             STOP

    END SELECT

    lda = n
    lwork = n
 
    allocate(ipiv(n))
    allocate(work(lwork))
    allocate(XT(n, size(vx(:,1))))
    allocate(X(size(vx(:,1)), n))
    allocate(XTX(n, n))
    
    iq = 0
    SELECT CASE (dimens)
      case (1)
! dimens = 1 
        do i = 0, d
            iq = iq + 1
            do ip = 1, size(vx(:,1))
               X(ip, i+1) = vx(ip, 1)**i  
            end do
        end do 

      case (2)
! dimens = 2 
        do i = 0, d
           do j = 0, d
                 if ( (i+j) .LE. d) then
                     iq = iq + 1
                     do ip = 1, size(vx(:,1))
                        X(ip, iq) = vx(ip, 1)**i    &
&                              * vx(ip, 2)**j    
                     end do
                 end if 
           end do
        end do

      case (3)
! dimens = 3 
        do i = 0, d
           do j = 0, d
              do k = 0, d
                 if ( (i+j+k) .LE. d) then
                     iq = iq + 1
                     do ip = 1, size(vx(:,1))
                        X(ip, iq) = vx(ip, 1)**i    &
&                                 * vx(ip, 2)**j    &
&                                 * vx(ip, 3)**k
                     end do
                 end if 
              end do
           end do
        end do
    END SELECT

!   ! prepare the matrix
!   do i = 0, d
!      do j = 1, size(vx)
!         X(j, i+1) = vx(j)**i
!      end do
!   end do
 
    XT  = transpose(X)
    XTX = matmul(XT, X)
 
    ! calls to LAPACK subs DGETRF and DGETRI
    call DGETRF(n, n, XTX, lda, ipiv, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if
    call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if
    polyfit = matmul( matmul(XTX, XT), vy)
 
    deallocate(ipiv)
    deallocate(work)
    deallocate(X)
    deallocate(XT)
    deallocate(XTX)
 
  end function
 
end module
