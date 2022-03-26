reset
#set terminal aqua
set terminal postscript eps enhanced 24
set output "fig4.eps"
set logscale xy
set xrange [10:1e-2]
set yrange [3e-14:3e-8]
set xlabel "T = temperature / MeV"
set ylabel "<{/Symbol s}v> = nuclear reaction rate / m^3 s^{-1}"
set xtics font "Helvetica,24"
set ytics font "Helvetica,24"
set format y "10^{%T}"
set key Left reverse at 0.1,3e-11 font "Helvetica,20"
plot 'fig4.txt' u 1:($2*1e6) t 'n + p -> d + {/Symbol g}' w l,\
     'fig4.txt' u 1:($4*1e6) t 'd + d -> ^3He + n' w l,\
     'fig4.txt' u 1:($5*1e6) t 'd + d -> t + p' w l,\
     'fig4.txt' u 1:($6*1e6) t '^3He + n -> t + p' w l,\
     'fig4.txt' u 1:($7*1e6) t 't + d -> ^4He + n' w l,\
     'fig4.txt' u 1:($8*1e6) t '^3He + d -> ^4He + p' w l
