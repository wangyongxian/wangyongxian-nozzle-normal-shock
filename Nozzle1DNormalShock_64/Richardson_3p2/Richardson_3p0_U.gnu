set grid
set time
set data style linespoints
set key right bottom
set xlabel 'log (h)'
set ylabel 'logaritmo do módulo de erros estimados'
  plot 'U_pL.dat'   using 1:2 title 'E(pL)'
replot 'U_p2.dat'   using 1:2 title 'E(p2)'
replot 'U_p3.dat'   using 1:2 title 'E(p3)'
replot 'U_12.dat'   using 1:2 title 'Ud(Th)'
replot 'U_12.dat'   using 1:3 title 'Uri_12(Th)'
replot 'Uri_pU.dat' using 1:2 title 'Uri_pU(Th)'
replot 'Uri_13.dat' using 1:2 title 'Uri_13(Th)'
replot 'Uri_23.dat' using 1:2 title 'Uri_23(Th)'
replot 'Uri_bi.dat' using 1:2 title 'Uri_bi(Th)'
replot 'Utri.dat'   using 1:2 title 'Uri_tri(Th)'
set title 'T(1/2), Peclet 1Dp 1.3, CDS, S(0), Pe=10, Richardson 3.1 (4/12/2007)'
replot
