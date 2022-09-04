reset session
set title "Trayectorias en el Eje Rotado" offset 0,-0.5
set xlabel "X_{rot}"
set ylabel  "Y_{rot}"

set title font "Helvetica, 18"
set xlabel font 'Helvetica,14'
set ylabel font 'Helvetica,14'
set key font "Helvetica, 11"

set xrange [-500:1100]
set yrange [-0.1 : 0.1]
set key top right box ls -1 lw 2
set grid
plot  "b_datos.txt" using 2:3 with lp ps 1 lw 2 t 'Sol', \
   "b_datos.txt" using 4:5 with lp ps 1 lw 2 t 'Jupiter' 
