reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig1.eps"
set logscale xy
set xrange [1e-2:1e4]
set yrange [1e-2:10]
set xlabel "t = time / sec"
set ylabel "T = temperature / MeV"
plot 'fig1.txt' u 1:2 t 'photon and electron' w l,\
     'fig1.txt' u 1:3 t 'neutrino' w l
