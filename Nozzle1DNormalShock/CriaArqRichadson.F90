Module CriaArqRichadson

implicit none
!modelo1.txt
!modelo2.txt
type ::ArqConf
character(90)   ::arq_dados
end type ArqConf

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


!subroutine 
!open(3,file="modelo1.txt")
!open(4,file="Richardson_3p0.in", )
!close(3,4)


!"D:\Projects\Richardson_3p0"
!"Richardson_3p0.exe"

!"Richardson_3p0.in" ! este arquivo diz qual eh o arquivo de configuracao para o programa
! exemplo da linha lida 'Richardson_3p0_Teste_001.in' 
!no arquivo de configuracao no caso:
!'Richardson_3p0_Teste_001.in' 
!linha:
!1 - o numero da variavel a ser analisada
!2 - nome do arquivo de saida
!3 - 0 sem solucao analitica 1 com solucao analitica
!4 - nome do arquivo da solucao analitica
!5 - ordem assintotica do erro verdadeiro
!6 - variação entre as ordens verdadeiras do erro verdadeiro
!7 - fator de segurança do estimador GCI
!8 -(tipo_refino) tipo de refino de malha (r):   1 = constante;   2 = variável
!9 -(numero_iteracoes) número máximo de iterações para obter a ordem aparente (pU) com r variável
!10 -(tol) tolerância na obtenção de pU (r variável)
!11 - (pU_min) limite inferior de pU (r variável)
!12 - (pU_max) limite superior de pU (r variável)   
!13 - (nm)  número de malhas
!14 - em seguida coloca-se o nomes dos arquivos com as solucoes de cada malha
!15 - e por ultimo o titulo ateh 90 caracteres

!arquivos de malha
!     4.0e-00      h (m)
! 1   2.048e+3     exato
! 2   2.560e+3     CDS
! 3   6.40e+2      DDS-2
! 4   1.152e+3     média
! 5   9.60e+2      UDS


!arquivo com a solucao analitica
! 1  2.048e+3   exato
! 2  2.048e+3   CDS
! 3  2.048e+3   DDS-2
! 4  8.192e+2   média
! 5  2.048e+3   UDS

subroutine WriteConfFile
!open(12, file="")

!close(12)
end subroutine WriteConfFile

subroutine CreateMeshFile

!open(12, file="")

!close(12)

end subroutine CreateMeshFile

end module CriaArqRichadson


