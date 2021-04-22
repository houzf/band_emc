set terminal postscript eps enhanced  color "Times-Roman" 25 
set output "bnd.eps"
set size ratio 0 0.668, 1
set border 31 lt -1 lw 3.000
set data style lines
set style line 1  lt 1 lc rgb "red"     lw 6.000 pt 1 ps default
set style line 11 lt 2 lc rgb "red"     lw 5.000 pt 3 ps default
set style line 2  lt 1 lc rgb "green"   lw 6.000 pt 1 ps default
set style line 22 lt 2 lc rgb "green"   lw 5.000 pt 3 ps default
set style line 3  lt 1 lc rgb "blue"    lw 6.000 pt 1 ps default
set style line 33 lt 2 lc rgb "blue"    lw 5.000 pt 3 ps default
set style line 4  lt 1 lc rgb "magenta" lw 6.000 pt 1 ps default
set style line 44 lt 2 lc rgb "magenta" lw 5.000 pt 3 ps default
set style line 5  lt 1 lc rgb "cyan"    lw 6.000 pt 1 ps default
set style line 55 lt 2 lc rgb "cyan"    lw 5.000 pt 3 ps default
set style line 8  lt 1 lc rgb "#696969" lw 3.000 pt 1 ps default
set style line 88 lt 2 lc rgb "#696969" lw 5.000 pt 3 ps default
set ticslevel 1.0 
set tics scale 2,1.0
set mytics 2
set ylabel "Energy (eV)" 0,0
set xtics ("G"  0.000000,  "X"  0.607988,  "X"  0.607988,  "W"  0.911983,  "W"  0.911983,  "L"  1.341895,  "L"  1.341895,  "G"  1.868429,  "G"  1.868429,  "K"  2.513298,  "K"  2.513298,  "X"  2.993955)
set xrange [  0.000000  :  2.993955]
set yrange [-12.000000  : 5.000000]
set arrow 3 from 0.0,  0.000000  to   2.993955  ,0.0 nohead lt 2 lw 4
plot "fatbnd-atom-1-spdf.dat" using 1:2 with lines ls 3 notitle,\
     "fatbnd-atom-1-spdf.dat" using 1:2:($3*10) with points lt 3 pt 3 ps variable notitle,\
     "fatbnd-atom-1-spdf.dat" using 1:2:($4*10) with points lt 4 pt 5 ps variable notitle,\
     "fatbnd-atom-1-spdf.dat" using 1:2:($5*10) with points lt 5 pt 7 ps variable notitle,\
     "highk.dat" using 1:2 with lines ls 8 notitle
