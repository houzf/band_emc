   SUBROUTINE write_highk(n, kdistances)
     use prec

     implicit none
     integer i, j
     integer, intent(in):: n
     real(dp), intent(in), dimension(n) :: kdistances 
     character(len=256) :: filename

     filename='highk.dat'
     open(10,file=filename)

     do i = 1, n
        do j = -8, 8 
          write(10,'(2F12.5)') kdistances(i), exp(0.5*j) &
 &         - exp(-0.5*j)
        end do
        write(10,*) 
     end do
     close(10)

     end subroutine write_highk
