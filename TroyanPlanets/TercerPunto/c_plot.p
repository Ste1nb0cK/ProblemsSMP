reset session
set title "Trayectorias en el Eje Rotado" offset 0,-0.5
set xlabel "X_{rot} [m_n]"
set ylabel  "Y_{rot} [m_n]"

set title font "Helvetica, 18"
set xlabel font 'Helvetica,14'
set ylabel font 'Helvetica,14'
set key font "Helvetica, 11"

set xrange [-200:1100]
set yrange [-20.1 : 990]
set key top right box ls -1 lw 2
set grid
plot  "c_datos.txt" using 2:3 with lp ps 2  lw 2 pt 22 t 'Sol', \
   "c_datos.txt" using 4:5 with lp ps 2 lw 2 pt 22 t 'Jupiter' , \
      "c_datos.txt" using 6:7 with lp ps 2 pt 22 lw 2 t ' Troyanos'
#set ytics offset -0.1,graph -0.06 rotate by 60
