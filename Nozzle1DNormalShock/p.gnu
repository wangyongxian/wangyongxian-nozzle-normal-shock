set data style linespoints
set grid
plot 'p.dat' using 2:3 title 'press�o'
set xlabel 'x (m)'
set ylabel 'p (Pa)'
set title 'PROG10_CFD1, M�TODO SIMPLEC, N=12 (5/10/2010)'                
replot
