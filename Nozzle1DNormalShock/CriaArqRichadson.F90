Module CriaArqRichadson

use variaveis

implicit none

contains

subroutine WriteConf1File(filename1, filename2)
    character*255, intent(in) ::filename1
    character*255, intent(in) ::filename2
character*255 form

open(14, file=filename1)
form = "'"//trim(adjustl(filename2)) //"'"
25 format(A52,'(arq_dados) nome do arquivo principal de dados do CASO a analisar (at� 50 caracteres)',1/)
write(14,25) form
close(14)
end subroutine WriteConf1File

subroutine WriteConf2File(filename1, filename2, filename3, filename_template, nmalhas)
    character*255, intent(in) ::filename1
    character*255, intent(in) ::filename2
    character*255, intent(in) ::filename3
    character*100, intent(in) ::filename_template
    integer     , intent(in)  ::nmalhas
    integer ::i
    character*255 ::form1
    character*255 ::form2
    character*255 ::arquivo
    
26 format (I8,  T30, '(variavel) n�mero da vari�vel, nos arquivos de resultados, a analisar', 1/, &
           A30, T30, '(caso)     nome do arquivo de sa�da', 1/, &
           I8,  T30, '(com_solucao_analitica) com solu��o anal�tica?   1 = sim;   0 = n�o', 1/, &
           A50, T50, '(arq_exato) nome do arquivo da solu��o anal�tica exata (at� 50 caracteres)', 1/, &
           I8,  T30, '(pL)  ordem assint�tica do erro verdadeiro', 1/, &
           I8,  T30, '(dpL) varia��o entre as ordens verdadeiras do erro verdadeiro', 1/, &
           d10.2,  T30, '(Fs)  fator de seguran�a do estimador GCI', 1/, &
           I8,  T30, '(tipo_refino) tipo de refino de malha (r):   1 = constante;   2 = vari�vel', 1/, &
           I8,  T30, '(numero_iteracoes) n�mero m�ximo de itera��es para obter a ordem aparente (pU) com r vari�vel', 1/, &
           d10.2,  T30, '(tol) toler�ncia na obten��o de pU (r vari�vel)', 1/, &
           d10.1,  T30, '(pU_min) limite inferior de pU (r vari�vel)', 1/, &
           d7.1,  T30, '(pU_max) limite superior de pU (r vari�vel)', 1/, &
           I8,  T30, '(nm)  n�mero de malhas')
open(12, file=filename1)
form1 = "'" // trim(adjustl(filename2)) // "'"
form2 = "'" // trim(adjustl(filename3)) // "'"
write(12,26) 1, form1 ,0, form2 ,2,2,&
                3.0d0, 1, 500,1.0d-15,-10.0d0,100.0d0,nmalhas
10 format(i8)
11 format(A50)

do i=1,nmalhas
write(arquivo,10) i
arquivo = trim(adjustl(filename_template)) // trim(adjustl(arquivo)) // ".Richardson_3p0"
files(i) = trim(adjustl(arquivo))
form1 = "'" // trim(adjustl(arquivo)) // "'"
form1 = trim(adjustl(form1))
write(12,11) form1
end do

12 format (A92, '(t�tulo) (at� 90 caracteres)',2/, &
                '012345    1         2         3         4         5         6         7         8         9        ')
title = "'" // trim(adjustl(title)) // "'"
write(12,12)  title  
close(12)

end subroutine WriteConf2File


subroutine WriteAnalitico(filename, throat, shock)
character*255, intent(in) ::filename
integer , intent(in) ::throat
real*8, intent(in) ::shock
integer ::i
i = throat
open(15,file=filename)
14 format (i2, T5 ,1PE25.16, T30 ,A15 )
!velocidade
write(15,14) 1, ue(i), 'velocidade - garganta'
!Mach
write(15,14) 2, (ue(i)/dsqrt(gama*R*0.5d0*(T(i)+T(i+1)))) , 'mach - garganta'
!perssao
write(15,14) 3, 0.5d0*(pa(i)+pa(i+1)), 'pressao - garganta'
!ro
write(15,14) 4, 0.5d0*(roa(i)+roa(i+1)), 'ro - garganta'
!temperatura
write(15,14) 5, 0.5d0*(ta(i)+ta(i+1)), 'temperatura - garganta'
!velocidade
write(15,14) 6, 0.5d0*(ua(1)+ua(2)), 'velocidade - entrada'
write(15,14) 7, 0.5d0*(ua(n-1)+ua(n)), 'velocidade - saida'
!pressao
write(15,14) 8, 0.5d0*(pa(1)+pa(2)), 'pressao - entrada'
write(15,14) 9, 0.5d0*(pa(n-1)+pa(n)), 'pressao - saida'
!temperatura
write(15,14) 10, 0.5d0*(ta(1)+ta(2)), 'temperatura - entrada'
write(15,14) 11, 0.5d0*(ta(n-1)+ta(n)), 'temperatura - saida'
!coeficiente de descarga
write(15,14) 12, 1.0d0, 'coeficiente de descarga - saida (adm)'
!empuxo dinamico
write(15,14) 13, 1.0d0, 'empuxo - saida (adm)'
write(15,14) 14, 0.0d0 , 'media da norma l1 - velocidade'
write(15,14) 15, 0.0d0 , 'media da norma l1 - temperatura'
write(15,14) 16, 0.0d0 , 'media da norma l1 - Ro'
write(15,14) 17, 0.0d0 , 'media da norma l1 - pressao'
write(15,14) 18, shock , 'local do choque (m)'

close(15)
end subroutine WriteAnalitico

subroutine WriteMeshFile(filename, throat, shock, h)
character*255, intent(in) ::filename
real*8, intent(in) ::h
integer, intent(in) ::throat
integer, intent(in) ::shock

open(16,file=filename)
12 format (T5 ,1pe25.16, T30 ,A15)
write(16,12) h, 'h (m)'
call WriteData(16,throat, shock)
close(16)


end subroutine WriteMeshFile

subroutine WriteData(descriptor, throat, shock)
integer, intent(in) ::descriptor
integer, intent(in) ::throat
integer, intent(in) ::shock
integer ::i
i = throat
14 format (i2, T5 ,1pe25.16, T50 ,A25)
!velocidade
write(descriptor,14) 1, ue(i), 'velocidade - garganta'
!Mach
write(descriptor,14) 2, (ue(i)/dsqrt(gama*R*0.5d0*(T(i)+T(i+1)))) , 'mach - garganta'
!perssao
write(descriptor,14) 3, 0.5d0*(p(i)+p(i+1)), 'pressao - garganta'
!ro
write(descriptor,14) 4, 0.5d0*(ro(i)+ro(i+1)), 'ro - garganta'
!temperatura
write(descriptor,14) 5, 0.5d0*(t(i)+t(i+1)), 'temperatura - garganta'
!velocidade
write(descriptor,14) 6, 0.5d0*(u(1)+u(2)), 'velocidade - entrada'
write(descriptor,14) 7, 0.5d0*(u(n-1)+u(n)), 'velocidade - saida'
!pressao
write(descriptor,14) 8, 0.5d0*(p(1)+p(2)), 'pressao - entrada'
write(descriptor,14) 9, 0.5d0*(p(n-1)+p(n)), 'pressao - saida'
!temperatura
write(descriptor,14) 10, 0.5d0*(t(1)+t(2)), 'temperatura - entrada'
write(descriptor,14) 11, 0.5d0*(t(n-1)+t(n)), 'temperatura - saida'
!coeficiente de descarga
write(descriptor,14) 12, roe(n-1)*ue(n-1)*Se(n-1)/(0.5d0*(Ma(n-1)+Ma(n))), 'coeficiente de descarga - saida (adm)'
!empuxo dinamico
write(descriptor,14) 13, roe(n-1)*ue(n-1)*Se(n-1)*ue(n-1)/(0.25d0*(Ma(n-1)+Ma(n))*(ua(n)+ua(n-1))), 'empuxo - saida (adm)'
write(descriptor,14) 14,sum(abs(Ua(2:n-1)-u(2:n-1)))/(N-2) , 'media da norma l1 - velocidade'
write(descriptor,14) 15,sum(abs(Ta(2:n-1)-T(2:n-1)))/(N-2) , 'media da norma l1 - temperatura'
write(descriptor,14) 16,sum(abs(Roa(2:n-1)-Ro(2:n-1)))/(N-2) , 'media da norma l1 - Ro'
write(descriptor,14) 17,sum(abs(Pa(2:n-1)-P(2:n-1)))/(N-2) , 'media da norma l1 - pressao'
if(P_OUT /= -1.0d0) then
write(descriptor,14) 18, (xp(shock+1)+xp(shock))/2.0d0 , 'local do choque'
end if

end subroutine WriteData


end module CriaArqRichadson