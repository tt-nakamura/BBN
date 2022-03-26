reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig8.eps"
set logscale x
set xrange [3e-11:1e-8]
set yrange [0.17:0.27]
set xlabel "{/Symbol h} = baryon to photon ratio"
set ylabel "X_4 = mass fraction of ^4He"
set format x "10^{%T}"
set key bottom
set label "T = 0.01MeV" at 2e-9,0.201
set label "{/Symbol t} = 885.7s" at 2e-9,0.195
plot 'fig8_n2.txt' u 1:7 t 'N_{/Symbol n} = 2' w l,\
     'fig7.txt'    u 1:7 t 'N_{/Symbol n} = 3' w l,\
     'fig8_n4.txt' u 1:7 t 'N_{/Symbol n} = 4' w l
