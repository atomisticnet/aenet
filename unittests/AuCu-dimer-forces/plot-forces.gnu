#!/usr/bin/env gnuplot

set encoding utf8

set style line  1 lt 1 lw 2.5 pt 6 ps 1.5 lc rgb "black"
set style line  2 lt 1 lw 2.5 pt 4 ps 1.5 lc rgb "red"
set style line  3 lt 1 lw 2.5 pt 8 ps 1.5 lc rgb "blue"
set style line  4 lt 1 lw 2.5 pt 6 ps 1.5 lc rgb "green"
set style line  5 lt 1 lw 2.5 pt 6 ps 1.5 lc rgb "orange"

set style line 11 lt 3 lw 3.5 pt 6 ps 1.5 lc rgb "black"
set style line 12 lt 3 lw 3.5 pt 4 ps 1.5 lc rgb "red"
set style line 13 lt 3 lw 3.5 pt 8 ps 1.5 lc rgb "blue"
set style line 14 lt 3 lw 3.5 pt 6 ps 1.5 lc rgb "green"
set style line 15 lt 3 lw 3.5 pt 6 ps 1.5 lc rgb "orange"

set style line 10 lt 1 lw 1.0 lc rgb "black"

set terminal postscript color eps enhanced size 5, 3 font "Times,24"

set border ls 10
set mxtics 2
set mytics 2

poly3(x,a0,a1,a2,a3) = a0 + a1*x + a2*x**2 + a3*x**3
dpoly3(x,a1,a2,a3) = a1+2.0*a2*x+3.0*a3*x**2

fit [*:3.0] poly3(x,a0,a1,a2,a3) 'forces.dat' u 1:2 via a0,a1,a2,a3

#----------------------------------------------------------------------#

set output 'graph-forces.eps'

set xlabel  "Atomic distance (Å)"
set ylabel  "Cohesive energy (eV)"
set y2label "Atomic force (eV/Å)"

set key bottom right

set ytics  nomirror format "%4.1f"
set y2tics nomirror format "%4.1f"

set yrange [-3.3:-2.8]
set y2range [-0.5:0.5]

set arrow 1 from 3.0, -3.3 to 3.0, -2.8 nohead ls 10

plot "forces.dat" u 1:2   w l ls  1 t "E(NN)", \
     poly3(x,a0,a1,a2,a3) w l ls 11 t "E(fit)", \
     "forces.dat" u 1:3   w l ls  2 t "F(NN)" axis x1y2, \
     -dpoly3(x,a1,a2,a3)  w l ls 12 t "F(deriv)" axis x1y2, \
     0 axis x1y2 w l ls 10 t ""

! epstopdf graph-forces.eps
! rm -f graph-forces.eps
