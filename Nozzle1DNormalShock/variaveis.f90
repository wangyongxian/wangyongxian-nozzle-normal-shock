module variaveis

use portlib

implicit none

integer :: N           ! n�mero total de n�s

integer :: i           ! n�mero do n�
                       ! i = 1, n� no contorno esquerdo fict�cio
                       ! i = N, n� no contorno direito fict�cio
                       ! 2 <= i <= N-1, n�s internos

real*8  :: tcpu        ! tempo de CPU em segundos

integer :: ver         ! aux�lio do comando System 

integer :: iteracao	   !	N�mero de itera��es

real*8  :: mi, ro      ! Viscosidade abs e massa espec�fica do fluido
real*8  :: Beta

real*8  :: fator, Uin  ! Fator de atrito e velocidade inicial

!real*8  :: Dzero, Cd   ! Di�metro inicial e coeficiente angular 
real*8  :: T0, Cd   ! Di�metro inicial e coeficiente angular 

real*8  :: Lt, deltax   ! Comprimento dom�nio de c�lculo e do volume de controle

real*8  :: deltat, Pi  ! N�mero de avan�os no tempo

real*8  :: Fat, R_o, R ! Coeficiente corre��o velocidade na face

real*8,dimension(:),allocatable :: u      ! solu��o num�rica
real*8,dimension(:),allocatable :: p      ! solu��o num�rica
real*8,dimension(:),allocatable :: rop      ! solu��o num�rica
real*8,dimension(:),allocatable :: roe      ! solu��o num�rica
real*8,dimension(:),allocatable :: plinha ! solu��o num�rica
real*8,dimension(:),allocatable :: Temperatura ! solu��o num�rica
real*8,dimension(:),allocatable :: u_o    ! solu��o num�rica inicial
real*8,dimension(:),allocatable :: ue     ! solu��o num�rica
real*8,dimension(:),allocatable :: ue_o   ! solu��o num�rica inicial
real*8,dimension(:),allocatable :: A, Ae  ! solu��o num�rica da �rea
real*8,dimension(:),allocatable :: M, Me  ! fluxo de massa na face leste

real*8,dimension(:),allocatable :: x, xe  ! coordenada espacial nodal
real*8,dimension(:),allocatable :: Raio      ! raio do duto

real*8 :: Rgases      ! constante dos gases
real*8 :: cp      ! cp
real*8 :: gama      ! gama
real*8 :: rin      ! rin
real*8 :: rg      ! rg
real*8 :: Lc      ! Lc
real*8 :: Ln      ! Ln


real*8,dimension(:),allocatable :: aPu, aPplinha ! coeficiente central de u e p
real*8,dimension(:),allocatable :: aWu, aWplinha ! coeficiente esquerdo de u e p
real*8,dimension(:),allocatable :: aEu, aEplinha ! coeficiente direito de u	e p
									   
real*8,dimension(:),allocatable :: bPu, bPplinha ! termo fonte de u e p

real*8,dimension(:),allocatable :: afu, atu, btu, bpru ! termos inclusos apu e bpu
real*8,dimension(:),allocatable :: ds, de   ! coeficientes do m�todo SIMPLEC

!character*30 :: nome      ! identifica��o do aluno
character*20 :: caso      ! nome do arquivo de sa�da

character*50 :: title     ! t�tulo do gr�fico
character*62 :: head      ! t�tulo do gr�fico + dia

character*12 :: dia       ! data da simula��o
character*8  :: hora      ! hor�rio da simula��o
integer*4    :: var(8)    ! data e hora
character*20 :: vardate   ! data e hora
character*20 :: vartime   ! data e hora
character*20 :: varzone   ! data e hora
character*70 :: note_caso ! notepad + caso
character*2  :: aux1,aux2
character*4  :: aux3
character*50 :: aux

end