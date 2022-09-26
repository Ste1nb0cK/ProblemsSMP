reset session
set title "Trayectorias Respecto al Centro de Masa" offset 0,-0.5
set xlabel "X [m_n]"
set ylabel  "Y [m_n]"

set title font "Helvetica, 18"
set xlabel font 'Helvetica,14'
set ylabel font 'Helvetica,14'
set key font "Helvetica, 11"

set xrange [-1200:1200]
set yrange [-1200:1200 ]
set key top right box ls -1 lw 2
set grid
plot  "a_datos.txt" using 1:2 with lp pt 22 ps 1 lw 2 t 'Sol', \
   "a_datos.txt" using 3:4 with l lw 2 t 'Jupiter' 
