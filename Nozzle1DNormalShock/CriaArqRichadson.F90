Module CriaArqRichadson

implicit none
!modelo1.txt
!modelo2.txt
type ::ArqConf
character(90)   ::arq_dados
end type ArqConf
character(120) ::linha
type :: ArqRichard
integer         ::variavel
character(90)   ::caso
integer         ::com_solucao_analitica
character(90)   ::arq_exato
integer         ::pl
integer         ::dpl
real*8          ::Fs
integer         ::tipo_refino
integer         ::numero_iteracoes
real*8          ::tol
real*8          ::pU_min
real*8          ::pU_max
integer         ::nm
character(90)   ::arq_numerico
character(90)   ::titulo
end type ArqRichard

contains

subroutine WriteConfFile(filename)
    character*50, intent(in) ::filename
open(12, file="modelo1.txt", MODE = 'READ')
read(12,*) linha
close(12)
open(12, file="Richardson_3p0.in")
write(12,*) '12', linha
close(12)
end subroutine WriteConfFile

subroutine CreateMeshFile(filename)
    character*50, intent(in) ::filename
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
           !'arq_numerico( 1) nome do arquivo com resultados numéricos na malha 1 (mais grossa)'
           !' arq_numerico( 2) nome do arquivo com resultados numéricos na malha 2'
open(12, file="teste.txt")
write(12,26) 5,"'Teste_030.txt'",0,"'Teste_analitico_001.Richardson_3p0'",1,1,&
                3.0d0, 1, 500,1.0d-15,-10.0d0,100.0d0,5
close(12)

end subroutine CreateMeshFile

end module CriaArqRichadson


!26 format (I8,  20x, A10, 2/, &
!           A10, 20x, A10, 2/, &
!           I8,  20x, A10, 2/, &
!           A50,  20x, A10, 2/, &
!           I8,  20x, A10, 2/, &
!           I8,  20x, A10, 2/, &
!           D14.3,  20x, A10, 2/, &
!           I8,  20x, A10, 2/, &
!           I8,  20x, A10, 2/, &
!           D14.3,  20x, A10, 2/, &
!           D14.3,  20x, A10, 2/, &
!           D14.3,  20x, A10, 2/, &
!           I8,  20x, A10, 2/)