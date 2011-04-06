set data style linespoints
set grid
plot 'fm.dat' using 2:3 title 'fluxo de massa'
set xlabel 'x (m)'
set ylabel 'fluxo de massa (kg/s)'
set title 'fluxo de massa'                
replot
