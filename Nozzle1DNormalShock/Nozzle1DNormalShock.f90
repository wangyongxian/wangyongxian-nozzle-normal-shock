!  Nozzle1DNormalShock.f90 
!
!  FUNCTIONS:
!  Nozzle1DNormalShock - Entry point of console application.
!

    program Nozzle1DNormalShock

!     Escoamento 1D de Fluido Incompress�vel 
!     Equa��o da MASSA e QML (Problema de Moody)


!     MODELO MATEM�TICO (resumo)
!       Equa��o diferencial: 
!          MASSA: d/dx(u*A) = 0
!          QML: ro*A*du/dt + ro*d/dx(A*(u**2)) = mi*d/dx(A*du/dx) - A*dp/dx - S(u,A) 
      
!       u = velocidade (inc�gnita)
!       p = press�o (inc�gnita) 
!       x = coordenada espacial 
!       S = termo fonte
!       Condi��es de contorno (C.C.) de Dirichlet

!     MODELO NUM�RICO (resumo)
!       M�todo num�rico: volumes finitos
!       malha uniforme unidimensional
!       solver: TDMA
!       precis�o: dupla
!       Linguagem Fortran
!       Aplicativo usado: Microsoft Developer Studio, Fortran PowerStation 4.0
!       Tipo de projeto: QuickWin Application
!       Tipos de aproxima��es num�ricas: 
!          u inc�gnita no termo advectivo: UDS
!          u coeficientes no termo difusivo: CDS
!       Acoplamento press�o-velocidade: M�todo SIMPLEC
!       Arranjo co-localizado de vari�veis
!       solu��o segregada das equa��es diferenciais
!       O processo � iterativo devido a n�o-linearidade da equa��o diferencial
!       os coeficientes dos contornos incorporam as C.C.
!          para u(0) = 0 (na entrada)
!          extrapola��o linear para u na sa�da
!       C.C. aplicadas com volumes fict�cios (P=1 e P=N)
!       propriedades constantes: ro(massa esp), f(fator de atrito), mi(viscosidade)
!       Express�o gen�rica do sistema de equa��es discretizado:
!          aP(i)*T(i) = aw(i)*T(i-1) + ae(i)*T(i+1) + bP(i) 
!              onde i = 1, 2, ... N (n�mero de n�s)

!     ARQUIVOS envolvidos no programa:
!       prog1_cdf1.f90 = programa principal
!       coef.f90       = calcula coeficientes e fontes do sistema linear
!       dados.f90      = l� e lista os dados do programa
!       result.f90     = resolve equa��es e gera listagens dos resultados
!       solvers.f90    = Solvers TDMA
!       variaveis.f90  = define todas as vari�veis globais do programa
!       prog5_cfd1.txt = arquivo de dados do programa
!       *.txt          = listagem dos resultados
!       T.dat          = arquivo de dados para fazer gr�fico
!       T.gnu          = arquivo de comandos para gerar gr�fico
!       notepad.exe    = editor dos arquivos
!       Wgnuplot.exe   = visualizador do gr�fico


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
	        //,'T�tulo = ', a<comp1>, &
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
