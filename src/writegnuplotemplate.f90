     subroutine gnuplot_template(n,labels,kdistances)
      use prec
      implicit none
      integer i 
      integer, intent(in):: n
      real(dp), intent(in), dimension(n) :: kdistances 
      character(len=1), intent(in), dimension(n) :: labels
      logical lexist
      character(len=256) :: filename

!.......... common part for both regimes ............................
      filename='test.gnu'
      inquire(file=filename,exist=lexist)
      if ( .not.  lexist) then
          write(*,'(3A)') "Error: ",trim(filename), " does not exist"
      else 
         call system('rm test.gnu')
      end if
      open (15, file='test.gnu', status='new')
      write (15,302)
 302  format ('set terminal postscript eps enhanced  color "Times-Roman" 25 '/ &
        'set output "bnd.eps"' /&
        'set size ratio 0 0.668, 1' /&
        'set border 31 lt -1 lw 3.000' /&
        'set data style lines'/ &
        'set style line 1  lt 1 lc rgb "red"     lw 6.000 pt 1 ps default'/&
        'set style line 11 lt 2 lc rgb "red"     lw 5.000 pt 3 ps default'/&
        'set style line 2  lt 1 lc rgb "green"   lw 6.000 pt 1 ps default'/&
        'set style line 22 lt 2 lc rgb "green"   lw 5.000 pt 3 ps default'/&
        'set style line 3  lt 1 lc rgb "blue"    lw 6.000 pt 1 ps default'/&
        'set style line 33 lt 2 lc rgb "blue"    lw 5.000 pt 3 ps default'/&
        'set style line 4  lt 1 lc rgb "magenta" lw 6.000 pt 1 ps default'/&
        'set style line 44 lt 2 lc rgb "magenta" lw 5.000 pt 3 ps default'/&
        'set style line 5  lt 1 lc rgb "cyan"    lw 6.000 pt 1 ps default'/&
        'set style line 55 lt 2 lc rgb "cyan"    lw 5.000 pt 3 ps default'/&
        'set style line 8  lt 1 lc rgb "#696969" lw 3.000 pt 1 ps default'/&
        'set style line 88 lt 2 lc rgb "#696969" lw 5.000 pt 3 ps default'/&
        'set ticslevel 1.0 '/ &
        'set tics scale 2,1.0'/ &
        'set mytics 2' / &
        'set ylabel "Energy (eV)" 0,0' &
             )
!    write(15,*)
    do i=1, n
      if (i>=n) then
         write(15,306, advance='no') labels(i), kdistances(i)
      else if (i>1) then
         write(15,305, advance='no') labels(i), kdistances(i) 
      else
         write(15,304, advance='no') labels(i), kdistances(i) 
      end if
    end do
 304 format('set xtics ("',A1,'"',F10.6,',',2X)
 305 format('"',A1,'"',F10.6,',',2X)
 306 format('"',A1,'"',F10.6,')')
     write(15,*)
     write(15,307) kdistances(1), kdistances(n)
 307 format('set xrange [',F10.6, 2X, ':', F10.6, ']') 
     write(15,308) -10.0, 10.0
 308 format('set yrange [',F10.6, 2X ':', F10.6, ']') 
     write(15,309) kdistances(1), kdistances(n) 
 309 format('set arrow 3 from 0.0,',F10.6, 2X, 'to ', F10.6, 2X, ',0.0 nohead lt 2 lw 4') 
     write (15,310)
 310 format('plot "bnd.dat" using 1:2 with lines ls 3 notitle,\'/ &
            '     "highk.dat" using 1:2 with lines ls 8 notitle' &
           )
    close (15)
!   call system ('gnuplot test.gnu')
end subroutine gnuplot_template 

