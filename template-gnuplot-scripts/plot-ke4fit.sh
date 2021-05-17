#!/bin/sh
if [ ! -f ke4fit.dat ] ; then
	echo "ke4fit.dat doesn't exist!"
	exit
fi

if [ -f inp.em ]; then
	nk=`cat inp.em | tail -1 | awk '{printf "%d\n", 2 * $2 + 1 }' `
	echo $nk
else
	echo "How many k points along each directions? Please enter an integer number"
	read nk
fi

#check command gsed
command -v  gsed >/dev/null 2>&1
if [ $? -eq 1 ]; then
	echo >&2 "gsed is not installed. Recommed to use gsed on MacOS. sed will be used"
	csed='sed'
else
	csed='gsed'
fi
#
for i  in `seq 1 7`
do
	case "$i" in
		1)
		dirt=100 ;;
		2)
		dirt=010 ;;
		3)
		dirt=001 ;;
		4)
		dirt=110 ;;
		5)
		dirt=101 ;;
		6)
		dirt=011 ;;
		7)
		dirt=111 ;;
		*)
		echo "unknown direction" ;;
	esac
fstl=`echo "scale=0; (($i-1) * $nk +  $i) " | bc `
lstl=`echo "scale=0; ($i * $nk + $i - 1) " | bc  `
echo $fstl   $lstl
sed -n "$fstl","$lstl"p  ke4fit.dat  > ke4fit-$dirt.dat
#rm -f fit.log
cat > plot-VB-ke-fit.plt <<!

reset

set encoding iso_8859_1
set terminal postscript eps enhanced color  "Times-Roman" 35
set output "VB-ke-$dirt-fit.eps"     # name of output file

# Color definitions
set border linewidth 3.5
set style line 1 lc rgb '#0060ad' lt 1 lw 2 # --- blue
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 7 # --- red

set size ratio 0 1.000, 1.3
# Axes
set style line 11 lc rgb '#808080' lt 1 lw 4
set border 3 back ls 11
set tics nomirror out scale 1.00
set arrow 11 from graph 1,0 to graph 1.05,0 size screen 0.025,15,60 filled ls 11
set arrow 22 from graph 0,1 to graph 0,1.05 size screen 0.025,15,60 filled ls 11
# Grid
set style line 12 lc rgb'#808080' lt -0 lw 1
set grid back ls 12

set style line 1  linetype 1 linecolor rgb "red"     linewidth 8.000 pointtype 1 pointsize default
set style line 11 linetype 2 linecolor rgb "red"     linewidth 5.000 pointtype 3 pointsize default
set style line 2  linetype 1 linecolor rgb "green"   linewidth 7.000 pointtype 1 pointsize default
set style line 22 linetype 2 linecolor rgb "green"   linewidth 5.000 pointtype 3 pointsize default
set style line 3  linetype 1 linecolor rgb "blue"    linewidth 8.000 pointtype 1 pointsize default
set style line 33 linetype 2 linecolor rgb "blue"    linewidth 5.000 pointtype 3 pointsize default
set style line 4  linetype 1 linecolor rgb "magenta"    linewidth 6.000 pointtype 1 pointsize default
set style line 44 linetype 2 linecolor rgb "magenta"    linewidth 5.000 pointtype 3 pointsize default
#set key at 5.2,1.5  noreverse enhanced box linetype -2 linewidth 1.200 samplen 2 spacing 1.5 width 0.0 height 0
set key right  top Right noreverse enhanced box linetype -2 linewidth 1.200 samplen 2 spacing 1.5 width 0.0 height 0  autotitles
set xlabel 'K-K_0 [1/Bohr]' offset 0,0.5
set ylabel 'Eigenvalue [Ha]' offset 1,0
set tics scale 1.50
set format y "%6.4f"
set mytics 5
set mxtics 5
set format x "%5.3f"
set xtics rotate by -45
#FIT_LIMIT= 1e-6
#set fit quiet
set fit logfile "$dirt-VB-fit.log"
f1(x)= a * x**2 + b * x + c
fit f1(x) "ke4fit-$dirt.dat" using 1:2 via a, b, c
f2(x)= p * x**2 + q * x + r
fit f2(x) "ke4fit-$dirt.dat" using 1:3 via p, q, r
f3(x)= u * x**2 + v * x + w
fit f3(x) "ke4fit-$dirt.dat" using 1:4 via u, v, w
#set datafile separator ","
plot f1(x) with lines ls 2  notitle '',\\
     "ke4fit-$dirt.dat" u 1:2   with points  lt 2 lw 4 pt 4 ps 4 notitle '' ,\\
     f2(x) with lines ls 3 notitle '', \\
     "ke4fit-$dirt.dat" u 1:3   with points  lt 3 lw 4 pt 6 ps 4 notitle '' ,\\
     f3(x) with lines ls 4 notitle '', \\
     "ke4fit-$dirt.dat" u 1:4   with points  lt 4 lw 4 pt 8 ps 4 notitle ''
!
/usr/bin/gnuplot  plot-VB-ke-fit.plt

a=`grep 'a      ' "$dirt"-VB-fit.log |grep '='| tail -1 | awk '{printf "%12.6f\n", $3}'`
p=`grep 'p      ' "$dirt"-VB-fit.log |grep '='| tail -1 | awk '{printf "%12.6f\n", $3}'`
u=`grep 'u      ' "$dirt"-VB-fit.log |grep '='| tail -1 | awk '{printf "%12.6f\n", $3}'`

#rm -f fit.log

cat > plot-CB-ke-fit.plt <<!

reset

set encoding iso_8859_1
set terminal postscript eps enhanced color  "Times-Roman" 35
set output "CB-ke-$dirt-fit.eps"     # name of output file

# Color definitions
set border linewidth 3.5
set style line 1 lc rgb '#0060ad' lt 1 lw 2 # --- blue
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 7 # --- red

set size ratio 0 1.000, 1.3
# Axes
set style line 11 lc rgb '#808080' lt 1 lw 4
set border 3 back ls 11
set tics nomirror out scale 1.00
set arrow 11 from graph 1,0 to graph 1.05,0 size screen 0.025,15,60 filled ls 11
set arrow 22 from graph 0,1 to graph 0,1.05 size screen 0.025,15,60 filled ls 11
# Grid
set style line 12 lc rgb'#808080' lt -0 lw 1
set grid back ls 12

set style line 1  linetype 1 linecolor rgb "red"     linewidth 8.000 pointtype 1 pointsize default
set style line 11 linetype 2 linecolor rgb "red"     linewidth 5.000 pointtype 3 pointsize default
set style line 2  linetype 1 linecolor rgb "green"   linewidth 7.000 pointtype 1 pointsize default
set style line 22 linetype 2 linecolor rgb "green"   linewidth 5.000 pointtype 3 pointsize default
set style line 3  linetype 1 linecolor rgb "blue"    linewidth 8.000 pointtype 1 pointsize default
set style line 33 linetype 2 linecolor rgb "blue"    linewidth 5.000 pointtype 3 pointsize default
set style line 4  linetype 1 linecolor rgb "magenta"    linewidth 6.000 pointtype 1 pointsize default
set style line 44 linetype 2 linecolor rgb "magenta"    linewidth 5.000 pointtype 3 pointsize default
#set key at 5.2,1.5  noreverse enhanced box linetype -2 linewidth 1.200 samplen 2 spacing 1.5 width 0.0 height 0
set key right  top Right noreverse enhanced box linetype -2 linewidth 1.200 samplen 2 spacing 1.5 width 0.0 height 0  autotitles
set xlabel 'K-K_0 [1/Bohr]' offset 0,0.5
set ylabel 'Eigenvalue [Ha]' offset 1,0
set tics scale 1.50
set format y "%6.4f"
set mytics 5
set mxtics 5
set format x "%5.3f"
set xtics rotate by -45
#FIT_LIMIT= 1e-6
#set fit quiet
set fit logfile "$dirt-CB-fit.log"
f1(x)= e * x**2 + f * x + g
fit f1(x) "ke4fit-$dirt.dat" using 1:5 via e, f, g
#set datafile separator ","
plot f1(x) with lines ls 4  notitle '',\\
     "ke4fit-$dirt.dat" u 1:5   with points  lt 4 lw 4 pt 4 ps 4 notitle ''
!
/usr/bin/gnuplot  plot-CB-ke-fit.plt

e=`grep 'e         ' "$dirt"-CB-fit.log | grep '=' |tail -1 | awk '{printf "%12.6f\n", $3}'`
if [ $i -eq  1  ]; then
  echo "$dirt" > EM-GNUPLOT.out
else
  echo "$dirt" >> EM-GNUPLOT.out
fi

if [  ${#a} -eq 0  ];then
a=10000
fi
if [  ${#p} -eq 0  ];then
p=10000
fi
if [  ${#u} -eq 0  ];then
u=10000
fi
if [  ${#e} -eq 0  ];then
e=10000
fi
echo $a, $p, $u, $e  | awk '{printf "%12.3f %12.3f %12.3f %12.3f\n", 0.5/$1, 0.5/$2, 0.5/$3, 0.5/$4}' >>EM-GNUPLOT.out
done
echo "# 0.000 --> 'Singular matrix in Invert_RtR' in the fitting of gnuplot." >>EM-GNUPLOT.out
echo "#Solution: manually fit it again using the command terminal of gnuplot instead of the script file." >>EM-GNUPLOT.out
