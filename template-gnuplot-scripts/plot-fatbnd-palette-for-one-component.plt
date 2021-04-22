# set this to the range of your variable which you want to color-encode
# or leave it out
set cbrange [0:1]
set yrange [-10:12]

# define the palette to your liking
set palette defined ( 0 "#B0B0B0", 0.333 "#FF0000", 0.666 "#0000FF", 1.0 "#000000" )

# in this example, column 3 is mapped to the colors of the palette
plot "fatbnd-tot-atoms-1-2.dat" u 1:2:3 w l lw 3  lc palette z
pause -1
