set grid
set time
set data style linespoints
set key right bottom
set xlabel 'log (h)'
set ylabel 'logaritmo do módulo de erros verdadeiros (E) e estimados (U)'
  plot 'E_pL_p2.dat' using 1:2 title 'E(pL)'
replot 'E_pL_p2.dat' using 1:3 title 'E(p2)'
replot 'E_p3.dat'   using 1:2 title 'E(p3)'
replot 'Eh.dat'     using 1:2 title 'E(Th)'
replot 'Ei_pU.dat'  using 1:2 title 'E(Ti_pU)'
replot 'Ebi.dat'    using 1:2 title 'E(Tbi_pU)'
replot 'U_12.dat'   using 1:3 title 'Uri_12(Th)'
replot 'U_12.dat'   using 1:4 title 'GCI_pL(Th)'
replot 'Uri_pU.dat' using 1:3 title 'Uri_p(Th)'
replot 'Uri_pU.dat' using 1:4 title 'GCI_p(Th)'
replot 'Uc_12.dat'  using 1:2 title 'Uc(Tc_12)'
replot 'Uc_13.dat'  using 1:2 title 'Uc(Tc_13)'
replot 'Ec_12.dat'  using 1:2 title 'E(Tc_12)'
replot 'Ec_13.dat'  using 1:2 title 'E(Tc_13)'
set title 'T(1/2), Peclet 1Dp 1.3, CDS, S(0), Pe=10, Richardson 3.2 (10/7/2009)'
replot
