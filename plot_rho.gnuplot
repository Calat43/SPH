#!/usr/bin/gnuplot 
set term gif \
    animate \
    optimize \
    delay 100 \
    size 512, 512 \
#    background "#ffeedf" \
#    crop \
#    font "Times-Roman,10"
set output "rho.gif"
set size square
set xrange[0:1]
set yrange[5:6]
set xzeroaxis
do for [i=0:499] {
	plot sprintf("rho/frame_%d.dat", i)
}
