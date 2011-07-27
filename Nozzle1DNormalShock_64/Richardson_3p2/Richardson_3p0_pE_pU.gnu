set grid
set time
set data style linespoints
set key right top
set yrange[-2:8]
set xlabel 'log (h)'
set ylabel 'ordens efetiva (pE) e aparente (pU)'
  plot 'pL.dat'     using 1:2 title 'pL'
replot 'p2.dat'     using 1:2 title 'pV(2)'
replot 'p3.dat'     using 1:2 title 'pV(3)'
replot 'pE_h.dat'   using 1:2 title 'pE(Th)'
replot 'pE_i.dat'   using 1:2 title 'pE(Ti_pU)'
replot 'pE_c12.dat' using 1:2 title 'pE(Tc_12)'
replot 'pE_c13.dat' using 1:2 title 'pE(Tc_13)'
replot 'pE_bi.dat'  using 1:2 title 'pE(Tbi_pU)'
replot 'pU_h.dat'   using 1:2 title 'pU(Th)'
replot 'pU_i.dat'   using 1:2 title 'pU(Ti_pU)'
replot 'pU_c12.dat' using 1:2 title 'pU(Tc_12)'
replot 'pU_c13.dat' using 1:2 title 'pU(Tc_13)'
replot 'pU_bi.dat'  using 1:2 title 'pU(Tbi_pU)'
set title 'T(1/2), Peclet 1Dp 1.3, CDS, S(0), Pe=10, Richardson 3.2 (10/7/2009)'
replot
