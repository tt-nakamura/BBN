reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig7.eps"
set logscale xy
set xrange [1e-11:1e-8]
set yrange [1e-6:1]
set xlabel "{/Symbol h} = baryon to photon ratio"
set ylabel "X = mass fraction"
set format xy "10^{%T}"
set key at 9e-9,2e-3
set label "T = 0.01MeV" at 1.5e-9,3e-2
set label "{/Symbol t} = 885.7s" at 1.5e-9,1.3e-2
set label "N_{/Symbol n} = 3" at 1.5e-9,5e-3
plot 'fig7.txt' u 1:3 t 'proton'  w l,\
     'fig7.txt' u 1:7 t '^4He' w l,\
     'fig7.txt' u 1:4 t 'deutron' w l,\
     'fig7.txt' u 1:6 t '^3He' w l,\
     'fig7.txt' u 1:5 t 'triton' w l
