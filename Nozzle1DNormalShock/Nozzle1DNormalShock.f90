!  Nozzle1DNormalShock.f90 
!
!  FUNCTIONS:
!  Nozzle1DNormalShock - Entry point of console application.
!

    program Nozzle1DNormalShock

!     Escoamento 1D de Fluido Incompressível 
!     Equação da MASSA e QML (Problema de Moody)


!     MODELO MATEMÁTICO (resumo)
!       Equação diferencial: 
!          MASSA: d/dx(u*A) = 0
!          QML: ro*A*du/dt + ro*d/dx(A*(u**2)) = mi*d/dx(A*du/dx) - A*dp/dx - S(u,A) 
      
!       u = velocidade (incógnita)
!       p = pressão (incógnita) 
!       x = coordenada espacial 
!       S = termo fonte
!       Condições de contorno (C.C.) de Dirichlet

!     MODELO NUMÉRICO (resumo)
!       Método numérico: volumes finitos
!       malha uniforme unidimensional
!       solver: TDMA
!       precisão: dupla
!       Linguagem Fortran
!       Aplicativo usado: Microsoft Developer Studio, Fortran PowerStation 4.0
!       Tipo de projeto: QuickWin Application
!       Tipos de aproximações numéricas: 
!          u incógnita no termo advectivo: UDS
!          u coeficientes no termo difusivo: CDS
!       Acoplamento pressão-velocidade: Método SIMPLEC
!       Arranjo co-localizado de variáveis
!       solução segregada das equações diferenciais
!       O processo é iterativo devido a não-linearidade da equação diferencial
!       os coeficientes dos contornos incorporam as C.C.
!          para u(0) = 0 (na entrada)
!          extrapolação linear para u na saída
!       C.C. aplicadas com volumes fictícios (P=1 e P=N)
!       propriedades constantes: ro(massa esp), f(fator de atrito), mi(viscosidade)
!       Expressão genérica do sistema de equações discretizado:
!          aP(i)*T(i) = aw(i)*T(i-1) + ae(i)*T(i+1) + bP(i) 
!              onde i = 1, 2, ... N (número de nós)

!     ARQUIVOS envolvidos no programa:
!       prog1_cdf1.f90 = programa principal
!       coef.f90       = calcula coeficientes e fontes do sistema linear
!       dados.f90      = lê e lista os dados do programa
!       result.f90     = resolve equações e gera listagens dos resultados
!       solvers.f90    = Solvers TDMA
!       variaveis.f90  = define todas as variáveis globais do programa
!       prog5_cfd1.txt = arquivo de dados do programa
!       *.txt          = listagem dos resultados
!       T.dat          = arquivo de dados para fazer gráfico
!       T.gnu          = arquivo de comandos para gerar gráfico
!       notepad.exe    = editor dos arquivos
!       Wgnuplot.exe   = visualizador do gráfico


! -----------------------------------------------

use dados

use resultados

! -----------------------------------------------

implicit none

integer :: comp, comp1

!-------------------------------------------------

  call date_and_time(vardate,vartime,varzone,var)

  write(aux,*) var(3)
  aux1 = trim(adjustl(aux))
  write(aux,*) var(2)
  aux2 = trim(adjustl(aux))
  write(aux,*) var(1)
  aux3 = trim(adjustl(aux))
  dia = '('//trim(aux1)//'/'//trim(aux2)//'/'//aux3//')'

  write(aux,*) var(5)
  aux1 = trim(adjustl(aux))
  write(aux,*) var(6)
  aux2 = trim(adjustl(aux))
  write(aux,*) var(7)
  aux3 = trim(adjustl(aux))
  hora = trim(aux1)//':'//trim(aux2)//':'//aux3

  call le_dados

  head = trim(title)//" "//trim(dia)//"'"

  open(10,file=caso)

  comp = len(trim(adjustl(nome)))

  comp1 = len(trim(adjustl(title)))

  write(10,18) trim(adjustl(nome)), trim(adjustl(title)), dia, hora
  18 format(/, 'Aluno  = ', a<comp>,  &
	        //,'Título = ', a<comp1>, &
            //,5x,'Dia = ',a12,5x,'Hora = ',a8)

  call mostra_dados 
  ! calcula a area e outras inicializacoes
  call inicializacao

  call solucao_numerica

  close (10)

  note_caso = 'notepad '//caso
  ver = system(note_caso) ! lista arquivo de resultados

! -----------------------------------------------

end
