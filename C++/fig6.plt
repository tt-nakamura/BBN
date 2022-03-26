reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig6.eps"
set logscale xy
set xrange [10:1e-2]
set yrange [0.08:1]
set xlabel "T = temperature / MeV"
set ylabel "X = mass fraction"
set key bottom left
set label "{/Symbol h} = 5x10^{-10}" at 7,0.21
set label "{/Symbol t} = 885.7s" at 7,0.175
set label "N_{/Symbol n} = 3" at 7,0.15
plot 'fig5-6.txt' u 1:2 t 'neutron' w l,\
     'fig5-6.txt' u 1:3 t 'proton'  w l,\
     'fig5-6.txt' u 1:7 t '^4He'  w l
