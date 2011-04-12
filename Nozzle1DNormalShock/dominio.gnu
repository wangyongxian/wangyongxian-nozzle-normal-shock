set data style linespoints
set grid
plot 'dominio.dat' using 1:2
set xlabel 'x (m)'
set ylabel 'y (m)'
set title 'dominio de calculo'                
replot
