#!/usr/bin/gnuplot 
set term gif \
    animate \
    optimize \
    delay 100 \
    size 512, 512 \
#    background "#ffeedf" \
#    crop \
#    font "Times-Roman,10"
set output "velocity.gif"
set size square
set xrange[0:1]
set yrange[-0.1:0.1]
set xzeroaxis
do for [i=0:499] {
	plot sprintf("velocity/frame_%d.dat", i)
}
