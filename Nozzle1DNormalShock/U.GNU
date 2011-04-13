set data style linespoints
set grid
plot 'u.dat' using 2:3 title 'velocidades' 
set xlabel 'x (m)'
set ylabel 'u (m/s)'
set title 'Resultado para N=12 (13/4/2011)'                              
replot
