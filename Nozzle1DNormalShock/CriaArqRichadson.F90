Module CriaArqRichadson

implicit none

contains

subroutine WriteConfFile(filename1, filename2)
    character*100, intent(in) ::filename1
    character*100, intent(in) ::filename2
character*100 form
open(14, file=filename1)
form = "'"//trim(adjustl(filename2)) //"'"
25 format(A50,'(arq_dados) nome do arquivo principal de dados do CASO a analisar (até 50 caracteres)',1/)
write(14,25) form
close(14)
end subroutine WriteConfFile

subroutine CreateMeshFile(filename1, filename2, filename3, filename_template, nmalhas)
    character*100, intent(in) ::filename1
    character*100, intent(in) ::filename2
    character*100, intent(in) ::filename3
    character*100, intent(in) ::filename_template
    integer     , intent(in)  ::nmalhas
    integer ::i
    character*100 ::form1
    character*100 ::form2
    character*100 ::arquivo
    
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
           I8,  T30, '(nm)  número de malhas', 1/)
open(12, file=filename1)
form1 = "'" // trim(adjustl(filename2)) // "'"
form2 = "'" // trim(adjustl(filename3)) // "'"
write(12,26) 5, form1 ,0, form2 ,1,1,&
                3.0d0, 1, 500,1.0d-15,-10.0d0,100.0d0,5
1000 format(i8)
do i=1, nmalhas
write(arquivo,1000) i
arquivo = trim(adjustl(arquivo)) // ".Richardson_3p0"
form1 = "'" // trim(adjustl(filename_template)) // arquivo
form1 = trim(adjustl(form1))
!'Teste_002.Richardson_3p0'  arq_numerico( 2) nome do arquivo com resultados numéricos na malha 2
write(12,*) form1
end do
close(12)

end subroutine CreateMeshFile

end module CriaArqRichadson

