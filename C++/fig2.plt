reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig2.eps"
set logscale xy
set xrange [1e2:1e-2]
set yrange [1e-10:1e10]
set xlabel "T = temperature / MeV"
set ylabel "{/Symbol G} = weak reaction rate / sec^{-1}"
set ytics 1000
set format y "10^{%T}"
plot 'fig2.txt' u 1:3 t 'neutron to proton' w l,\
     'fig2.txt' u 1:2 t 'proton to neutron' w l,\
     'fig2.txt' u 1:4 t 'cosmic expansion rate' w l
