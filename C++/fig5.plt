reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig5.eps"
set logscale xy
set xrange [2:1e-2]
set yrange [1e-24:2]
set xlabel "T = temperature / MeV"
set ylabel "X = mass fraction"
set ytics 1e5
set format y "10^{%T}"
set key bottom
set label "{/Symbol h} = 5x10^{-10}" at 1.7,2e-15
set label "{/Symbol t} = 885.7s" at 1.7,4e-17
set label "N_{/Symbol n} = 3" at 1.7,1e-18
plot 'fig5-6.txt' u 1:2 t 'neutron' w l,\
     'fig5-6.txt' u 1:3 t 'proton'  w l,\
     'fig5-6.txt' u 1:4 t 'deutron' w l,\
     'fig5-6.txt' u 1:5 t 'triton' w l,\
     'fig5-6.txt' u 1:6 t ' ^3He' w l,\
     'fig5-6.txt' u 1:7 t ' ^4He' w l,\
     'fig5-6.txt' u 1:8 t ' ^7Li' w l,\
     'fig5-6.txt' u 1:9 t ' ^7Be' w l
