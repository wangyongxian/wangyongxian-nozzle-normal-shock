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
25 format(A52,'(arq_dados) nome do arquivo principal de dados do CASO a analisar (até 50 caracteres)',1/)
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
    
26 format (I8,  T30, '(variavel) número da variável, nos arquivos de resultados, a analisar', 1/, &
           A30, T30, '(caso)     nome do arquivo de saída', 1/, &
           I8,  T30, '(com_solucao_analitica) com solução analítica?   1 = sim;   0 = não', 1/, &
           A50, T50, '(arq_exato) nome do arquivo da solução analítica exata (até 50 caracteres)', 1/, &
           I8,  T30, '(pL)  ordem assintótica do erro verdadeiro', 1/, &
           I8,  T30, '(dpL) variação entre as ordens verdadeiras do erro verdadeiro', 1/, &
           d10.2,  T30, '(Fs)  fator de segurança do estimador GCI', 1/, &
           I8,  T30, '(tipo_refino) tipo de refino de malha (r):   1 = constante;   2 = variável', 1/, &
           I8,  T30, '(numero_iteracoes) número máximo de iterações para obter a ordem aparente (pU) com r variável', 1/, &
           d10.2,  T30, '(tol) tolerância na obtenção de pU (r variável)', 1/, &
           d10.1,  T30, '(pU_min) limite inferior de pU (r variável)', 1/, &
           d7.1,  T30, '(pU_max) limite superior de pU (r variável)', 1/, &
           I8,  T30, '(nm)  número de malhas')
open(12, file=filename1)
form1 = "'" // trim(adjustl(filename2)) // "'"
form2 = "'" // trim(adjustl(filename3)) // "'"
write(12,26) 5, form1 ,0, form2 ,1,1,&
                3.0d0, 1, 500,1.0d-15,-10.0d0,100.0d0,5
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

12 format (A92, '(título) (até 90 caracteres)',2/, &
                '012345    1         2         3         4         5         6         7         8         9        ')
title = "'" // trim(adjustl(title)) // "'"
write(12,12)  title  
close(12)

end subroutine WriteConf2File


subroutine WriteAnalitico(filename, throat)
character*255, intent(in) ::filename
integer , intent(in) ::throat
integer ::i
i = throat
open(15,file=filename)
14 format (i2, T5 ,1PE21.13, T30 ,A15 )
!velocidade
write(15,14) 1, ue(i), 'velocidade'
!Mach
write(15,14) 2, (ua(i)/dsqrt(gama*R*T(i))) , 'mach'
!perssao
write(15,14) 3, pa(i), 'pressao'
!ro
write(15,14) 4, roa(i), 'roa'
!temperatura
write(15,14) 5, ta(i), 'temperatura'
close(15)
end subroutine WriteAnalitico

subroutine WriteMeshFile(filename, throat, h)
character*255, intent(in) ::filename
real*8, intent(in) ::h
integer, intent(in) ::throat

open(16,file=filename)
12 format (T5 ,1pe21.13, T30 ,A15)
write(16,12) h, 'h (m)'
call WriteData(16,throat)
close(16)


end subroutine WriteMeshFile

subroutine WriteData(descriptor, throat)
integer, intent(in) ::descriptor
integer, intent(in) ::throat
integer ::i
i = throat
14 format (i2, T5 ,1pe21.13, T30 ,A15)
!velocidade
write(descriptor,14) 1, ue(i), 'velocidade'
!Mach
write(descriptor,14) 2, (u(i)/dsqrt(gama*R*T(i))) , 'mach'
!perssao
write(descriptor,14) 3, p(i), 'pressao'
!ro
write(descriptor,14) 4, ro(i), 'ro'
!temperatura
write(descriptor,14) 5, t(i), 'temperatura'

end subroutine WriteData

!USE IFPORT
!result = CHDIR(dir_name)


end module CriaArqRichadson