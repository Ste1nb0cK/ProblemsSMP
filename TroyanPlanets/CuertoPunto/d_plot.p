reset session
set title "Trayectoria Perturbada del Troyano" offset 0,-0.5
set xlabel "X_{rot} [m_n]"
set ylabel  "Y_{rot} [m_n]"

set title font "Helvetica, 18"
set xlabel font 'Helvetica,14'
set ylabel font 'Helvetica,14'
set key font "Helvetica, 11"

set xrange [460:540]
set yrange [840 : 890]
set key top right box ls -1 lw 2
set grid
plot "d_datos.txt" using 6:7 with l lw 2 t ' Troyano'
#set ytics offset -0.1,graph -0.06 rotate by 60
