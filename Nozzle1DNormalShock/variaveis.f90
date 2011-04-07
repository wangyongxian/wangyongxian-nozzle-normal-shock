module variaveis

use portlib

implicit none

integer :: N           ! número total de nós

integer :: i           ! número do nó
                       ! i = 1, nó no contorno esquerdo fictício
                       ! i = N, nó no contorno direito fictício
                       ! 2 <= i <= N-1, nós internos

real*8  :: tcpu        ! tempo de CPU em segundos

integer :: ver         ! auxílio do comando System 

integer :: iteracao	   !	Número de iterações

real*8  :: mi, ro      ! Viscosidade abs e massa específica do fluido
real*8  :: Beta

real*8  :: fator, Uin  ! Fator de atrito e velocidade inicial

!real*8  :: Dzero, Cd   ! Diâmetro inicial e coeficiente angular 
real*8  :: T0, Cd   ! Diâmetro inicial e coeficiente angular 

real*8  :: Lt, deltax   ! Comprimento domínio de cálculo e do volume de controle

real*8  :: deltat, Pi  ! Número de avanços no tempo

real*8  :: Fat, R_o, R ! Coeficiente correção velocidade na face

real*8,dimension(:),allocatable :: u      ! solução numérica
real*8,dimension(:),allocatable :: p      ! solução numérica
real*8,dimension(:),allocatable :: rop      ! solução numérica
real*8,dimension(:),allocatable :: roe      ! solução numérica
real*8,dimension(:),allocatable :: plinha ! solução numérica
real*8,dimension(:),allocatable :: Temperatura ! solução numérica
real*8,dimension(:),allocatable :: u_o    ! solução numérica inicial
real*8,dimension(:),allocatable :: ue     ! solução numérica
real*8,dimension(:),allocatable :: ue_o   ! solução numérica inicial
real*8,dimension(:),allocatable :: A, Ae  ! solução numérica da Área
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
real*8,dimension(:),allocatable :: ds, de   ! coeficientes do método SIMPLEC

!character*30 :: nome      ! identificação do aluno
character*20 :: caso      ! nome do arquivo de saída

character*50 :: title     ! título do gráfico
character*62 :: head      ! título do gráfico + dia

character*12 :: dia       ! data da simulação
character*8  :: hora      ! horário da simulação
integer*4    :: var(8)    ! data e hora
character*20 :: vardate   ! data e hora
character*20 :: vartime   ! data e hora
character*20 :: varzone   ! data e hora
character*70 :: note_caso ! notepad + caso
character*2  :: aux1,aux2
character*4  :: aux3
character*50 :: aux

end