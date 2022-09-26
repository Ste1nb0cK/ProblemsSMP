reset session
set title "Componente X del Troyano Perturbado" offset 0,-0.5
set ylabel "X_{rot} [m_n]"
set xlabel  "Tiempo [s_n]"

set title font "Helvetica, 18"
set xlabel font 'Helvetica,14'
set ylabel font 'Helvetica,14'
set key font "Helvetica, 11"

#set arrow from 39806, graph 0 to 39806, graph 1 nohead
set xrange [0 : 78000]
set yrange [460 : 540]
set key top right box ls -1 lw 2
set grid
plot "d_datos.txt" using 1:6 with l lw 2 t ' Troyano'
#set ytics offset -0.1,graph -0.06 rotate by 60
