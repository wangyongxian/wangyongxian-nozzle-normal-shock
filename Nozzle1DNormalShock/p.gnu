set data style linespoints
set grid
plot 'p.dat' using 2:3 title 'press�o'
set xlabel 'x (m)'
set ylabel 'p (Pa)'
set title 'pressao'                
replot
